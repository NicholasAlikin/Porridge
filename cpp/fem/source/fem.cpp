#include "fem.hpp"

namespace fem {
size_t Point::EmptyID = 0;
size_t BaseElement::EmptyID = 0;

bool BaseNode::operator<(const BaseNode& other) const {
	return id < other.id;
}

GlobalDofs::GlobalDofs(size_t ndofs, const math::vector_t<size_t,2>& ID)
		: ndofs(ndofs), ID(ID) {};
GlobalDofs::GlobalDofs(size_t ndofs, math::vector_t<size_t,2>&& ID)
		: ndofs(ndofs), ID(std::move(ID)) {};

GlobalDofs::GlobalDofs(const GlobalDofs& other)
		: GlobalDofs(other.ndofs,other.ID) {};
GlobalDofs::GlobalDofs(GlobalDofs&& other)
		: GlobalDofs(other.ndofs, std::move(other.ID)) {};

BaseElement::BaseElement(const BaseElement& other)
		: ID(other.ID),propID(other.propID),matlID(other.matlID)
		 ,nodes(other.nodes),parameters(other.parameters) {}
BaseElement::BaseElement(BaseElement&& other)
		: ID(other.ID),propID(other.propID),matlID(other.matlID)
		 ,nodes(std::move(other.nodes)),parameters(std::move(other.parameters)) {}
BaseElement::BaseElement(size_t ID, size_t propID, size_t matlID
					   , const math::vector<size_t> &nodes
					   , const math::vector<double> &parameters)
		: ID(ID),propID(propID),matlID(matlID)
		 ,nodes(nodes),parameters(parameters) {}

size_t BaseElement::ndofs() const {
	return nnodes()*ndofs_node();
}



void BaseElement::set_nodes(const math::vector<size_t>& nodes) {
	if (nodes.size() != nnodes())
		std::logic_error("Incorrect number of element nodes!");
	this->nodes = nodes;
}
void BaseElement::set_nodes(math::vector<size_t>&& nodes) {
	if (nodes.size() != nnodes()) {	
		std::logic_error("Incorrect number of element nodes!");
	}
	this->nodes = std::move(nodes);
}

GlobalDofs parse_nodal_data(const math::vector<Node>& nodes_info
                          , const math::vector<NodeConstraint>& constraints_info
						  , math::vector<BaseElement*>& elements)
{
    // matrix with nodes global dofs
    math::vector_t<size_t,2> ID = math::zeros<size_t>(nodes_info.size(),BaseNode::DOFS);
    
    // iterators
    auto gnode = ID.begin(), node_end = ID.end();
    typename math::vector<size_t>::iterator gdof, gdof_end;

    auto node_info = nodes_info.begin();
    auto constraint_node = constraints_info.begin();
	typename decltype(NodeConstraint::dofs)::const_iterator constraint_dof;
    
	// set ID matrix by global dofs if current dof is not constrained
	// constrained dofs value = 0
    size_t global_id = 1;
    while (gnode != node_end) {
        gdof = gnode->begin();
		gdof_end = gnode->end();
        
		if (node_info->id == constraint_node->id) {
			constraint_dof = constraint_node->dofs.begin();
			while (gdof != gdof_end) {
				if (*constraint_dof == 0) {
					*gdof = global_id;
					++global_id;
				}
				++gdof; ++constraint_dof;
			}
			++constraint_node;
        } else {
			while (gdof != gdof_end) {
				*gdof = global_id;
				++global_id; ++gdof;
			}
        }
		++gnode; ++node_info;
    }
	
	// Calculate elements properties
	math::vector<Node> elem_nodes;
	decltype(elem_nodes.begin()) node;
	
	auto elem = elements.begin(), elem_end = elements.end();
	while (elem != elem_end) {
		elem_nodes = math::zeros<Node>((*elem)->nnodes());
		node = elem_nodes.begin();
		for (size_t nodeID: (*elem)->nodes) {
			*node = nodes_info[nodeID];
			++node;
		}
		(*elem)->calc_parameters(elem_nodes);
		++elem;
	}
	
	return {global_id-1, std::move(ID)};
}

void store_load_vector(math::vector<double> load, const GlobalDofs& GDofs
									 , const math::vector<NodeLoad>& loads_info)
{
	auto load_node     = loads_info.begin()
		,load_node_end = loads_info.end();
		
	decltype(loads_info[0].dofs.begin()) load_dof, load_dof_end;
	
	decltype(GDofs.ID.begin()) gnode;
	decltype(GDofs.ID[0].begin()) gdof;
	
	while (load_node != load_node_end) {
		load_dof = load_node->dofs.begin();
		load_dof_end = load_node->dofs.end();
		gdof = GDofs.ID[load_node->id].begin();
		while (load_dof != load_dof_end) {
			if (*gdof > 0) {
				load[*gdof-1] = *load_dof;
			}
			++load_dof; ++gdof;
		}
		++load_node;
	}
}



void assemble_precomputing(GlobalDofs& GDofs, const math::vector<BaseElement*>& elements) {
	/*
	Assemble global matrices by two loops
	1. define elements global dofs and global matrix column heights
	2. fill global matrix
	
	Store global matrix as its upper triangular part
	*/
	
	// Elements global dofs
	math::vector_t<size_t,2> ElemsGDofs(elements.size()); // might be not rectangular matrix!
	// Global matrix columns heights
	math::vector<size_t> colhs(GDofs.ndofs);
	size_t min_dof;
	
	// 1-st Loop over all elements
	auto elem = elements.begin(), elem_end = elements.end();
	auto elemdofs = ElemsGDofs.begin(); // iterator oncurrent element dofs
	decltype(ElemsGDofs[0].begin()) elemdof, elemdof_end; // iterator on current element dof
	decltype(elements[0]->nodes.begin()) elemnode, elemnode_end; // iterator on element node
	decltype(GDofs.ID[0].begin()) nodegdof, nodegdof_end; // iterator on element global dof
	while (elem != elem_end) {
		*elemdofs = math::zeros<size_t>((*elem)->ndofs()); // init elemdofs array
		elemdof = elemdofs->begin();
		elemnode = (*elem)->nodes.begin();
		elemnode_end = (*elem)->nodes.end();
		// loop over element nodes
		while (elemnode != elemnode_end) {
			nodegdof = GDofs.ID[*elemnode].begin();
			nodegdof_end = nodegdof + (*elem)->ndofs_node();
			while (nodegdof != nodegdof_end) {
				*elemdof = *nodegdof;
				++elemdof;
				++nodegdof;
			}

			++elemnode;
		}
		// calc current column heights
		min_dof = fem::min_elem_dof(elemdofs->begin(),elemdofs->end());
		elemdof = elemdofs->begin();
		elemdof_end = elemdofs->end();
		while (elemdof != elemdof_end) {
			if (*elemdof+1 > min_dof) {
				colhs[*elemdof-1] = std::max(*elemdof+1 - min_dof,colhs[*elemdof-1]);
			}
			++elemdof;
		}

		++elem; ++elemdofs;
	}


	math::vector<size_t> diags(GDofs.ndofs+1);
	auto diag = diags.begin(), diag_end = diags.end()-1;
	auto colh_1 = colhs.begin();
	*diag = 0;
	++diag;
	while (diag != diag_end) {
		*diag = *(diag-1) + *colh_1;
		++diag;	++colh_1;
	}
	diags.last() = diags[GDofs.ndofs-1] + *colh_1;


	GDofs.stiffness = math::zeros<double>(sum(colhs));
	GDofs.ElemID = std::move(ElemsGDofs);
	GDofs.diags = std::move(diags);
	// size_t band_width = *std::max_element(colhs.begin(),colhs.end());
	GDofs.colhs = std::move(colhs);


}

math::vector_t<double,2> assemble(GlobalDofs& GDofs, const math::vector<Elem*>& elements
				, const math::vector_t<double,2>& properties
				, const math::vector_t<double,2>& materials
				, math::vector_t<double>& band
				, math::vector_t<double,2> (Elem::* element_matrix)(const math::vector<double>&,const math::vector<double>&) const)
{
	// 2. fill global matrix
	// TODO
	math::vector_t<double,2> matrix = math::zeros<double>(GDofs.ndofs,GDofs.ndofs);
	// band = math::zeros<double>(math::sum(GDofs.colhs));
	
	math::vector_t<double,2> matrix_local; // element local matrix
	decltype(matrix_local.begin()) matloc_row,matloc_row_end;
	decltype(matrix_local[0].begin()) matloc_col,matloc_col_end;

	auto elem = elements.begin(), elem_end = elements.end();
	auto elemgdofs = GDofs.ElemID.begin(); // iterator on current element dofs
	decltype(GDofs.ElemID[0].begin()) elemgdof_row, elemgdof_row_end
									, elemgdof_col, elemgdof_col_end; // iterator on current element dof
	

	decltype(elements[0]->nodes.begin()) elemnode, elemnode_end; // iterator on element node
	decltype(GDofs.ID[0].begin()) nodegdof, nodegdof_end; // iterator on element global dof
	size_t matloc_row_num;
	// loop over elements
	while (elem != elem_end) {
		// calculate element matrix (stiffness, mass)
		matrix_local = ((*elem)->*element_matrix)(properties[(*elem)->propID]
										         , materials[(*elem)->matlID]);
		
		// loop over element matrix
		// loop over matrix rows
		for (matloc_row = matrix_local.begin()    // iterator on elem matrix row
		    ,matloc_row_end = matrix_local.end()  // same
		    ,elemgdof_row = elemgdofs->begin()    // iterator on global dofs, corresponding row global dofs of the elem
			,matloc_row_num = 0					  // elem matrix count - to go throw the upper triangular part of the matrix only
		   				;matloc_row != matloc_row_end // loop over all matrix rows
									;++matloc_row
									,++elemgdof_row
									,++matloc_row_num) {
			if (*elemgdof_row == 0) continue;     // assamble only matrix rows, which correspond not constrained global dofs

			// loop over columnt of current row
			for (matloc_col = matloc_row->begin() + matloc_row_num  // iterator on column element - only upper triangular part
			    ,matloc_col_end = matloc_row->end()                 // same
				,elemgdof_col = elemgdofs->begin() + matloc_row_num // iterator on global dofs, corresponding col global dofs of the elem
							;matloc_col != matloc_col_end // loop over row elements from diagonal to the end
										;++matloc_col
										,++elemgdof_col) {
				if (*elemgdof_col == 0) continue;  // assamble only matrix columns, which correspond not constrained global dofs


				// global matrix like 2d array - stored only upper triangular part, lower triangular is symmetric
				matrix[*elemgdof_row-1][*elemgdof_col-1] += *matloc_col;
				
				// global matrix like 1d array - stored only upper triangular part
				// GDofs.diags - addresses of the matrix diagonal elements
				// &Kij = &diag + col_num - row_num
				band[(GDofs.diags[*elemgdof_col-1]
				     + *elemgdof_col) - *elemgdof_row] += *matloc_col;
				
			}
		}
		
		
		++elem; ++elemgdofs;
	}

	return matrix;
}






} // namespace fem