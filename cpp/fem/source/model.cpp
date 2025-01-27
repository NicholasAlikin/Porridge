#include "model.hpp"

namespace fem {

/*==========
    Model
===========*/
Model::Model(const Model& other)
        :nodes_info(other.nodes_info)
        ,loads_info(other.loads_info)
        ,loads_var_info(other.loads_var_info)
		,constraints_info(other.constraints_info)
        ,elements(other.elements)
        ,materials(other.materials)
        ,properties(other.properties)
        ,elements_releases(other.elements_releases)
        ,nodes_releases(other.nodes_releases)
{}
Model::Model(Model&& other)
        :nodes_info(std::move(other.nodes_info))
        ,loads_info(std::move(other.loads_info))
        ,loads_var_info(std::move(other.loads_var_info))
		,constraints_info(std::move(other.constraints_info))
        ,elements(std::move(other.elements))
        ,materials(std::move(other.materials))
        ,properties(std::move(other.properties))
        ,elements_releases(std::move(other.elements_releases))
        ,nodes_releases(std::move(other.nodes_releases))
{}

Model::Model(const math::vector<Node>& nodes_info					
	 		,const math::vector<NodeLoad>& loads_info				
	 		,const math::vector<NodeConstraint>& constraints_info	
	 		,const math::vector<BaseElement*>& elements			
	 		,const math::vector_t<double,2>& materials				
	 		,const math::vector_t<double,2>& properties			
	 		,const math::vector<ElemReleases>& elements_releases)
        :nodes_info(nodes_info)
        ,loads_info(loads_info)
        ,constraints_info(constraints_info)
        ,elements(elements)
        ,materials(materials)
        ,properties(properties)
        ,elements_releases(elements_releases)
        ,nodes_releases(releases_precomputing(Model::elements_releases,elements))
{}

Model::Model(const math::vector<Node>& nodes_info					
	 		,const math::vector<NodeLoadVar>& loads_info				
	 		,const math::vector<NodeConstraint>& constraints_info	
	 		,const math::vector<BaseElement*>& elements			
	 		,const math::vector_t<double,2>& materials				
	 		,const math::vector_t<double,2>& properties			
	 		,const math::vector<ElemReleases>& elements_releases)
        :nodes_info(nodes_info)
        ,loads_var_info(loads_info)
        ,constraints_info(constraints_info)
        ,elements(elements)
        ,materials(materials)
        ,properties(properties)
        ,elements_releases(elements_releases)
        ,nodes_releases(releases_precomputing(Model::elements_releases,elements))
{}

Assemble::Assemble(const Assemble& other)
		:ndofs(other.ndofs)
		,nodes_dofs(other.nodes_dofs)
		,elems_dofs(other.elems_dofs)
		,colhs(other.colhs)
		,diags(other.diags)
		,band_width(other.band_width)
		,max_elem_dofs(other.max_elem_dofs)
		,elem_matrix(other.elem_matrix)
		,elem_load(other.elem_load)
		,elem_displ(other.elem_displ)
{}

Assemble::Assemble(Assemble&& other)
		:ndofs(other.ndofs)
		,nodes_dofs(std::move(other.nodes_dofs))
		,elems_dofs(std::move(other.elems_dofs))
		,colhs(std::move(other.colhs))
		,diags(std::move(other.diags))
		,band_width(other.band_width)
		,max_elem_dofs(std::move(other.max_elem_dofs))
		,elem_matrix(std::move(other.elem_matrix))
		,elem_load(std::move(other.elem_load))
		,elem_displ(std::move(other.elem_displ))
{}

/* Count all releases in nodes */
math::vector<NodeAllReleases> Model::releases_precomputing(math::vector<ElemReleases>& elements_releases
											       , const math::vector<BaseElement*>& elements) {
	
	
	std::unordered_map<size_t,NodeAllReleases> node_releases;
	auto elem_release = elements_releases.begin();
	auto elem_release_end = elements_releases.end();
	decltype(elements_releases[0].nodes.begin()) node_release, node_release_end;
	decltype(NodeAllReleases::dofs)* gnode_dofs;
	decltype(gnode_dofs->begin()) gnode_dof;
	decltype(elements_releases[0].nodes[0].dofs.begin()) node_dof, node_dof_end;

	size_t node_global_id;
	// loop over all elements with released nodes
	while (elem_release != elem_release_end) {
		node_release = elem_release->nodes.begin();
		node_release_end = elem_release->nodes.end();
		
		
		// loop over all element released nodes
		while (node_release != node_release_end) {
			node_global_id = elements[elem_release->elem_id]->nodes[node_release->node_localid];
			gnode_dofs = &node_releases[node_global_id].dofs;

			node_release->shift = node_releases[node_global_id].dofs.size()/BaseNode::DOFS;

			gnode_dofs->resize((node_release->shift+1)*BaseNode::DOFS);

			gnode_dof = gnode_dofs->begin() + node_release->shift*BaseNode::DOFS;
			node_dof = node_release->dofs.begin();
			node_dof_end = node_release->dofs.end();

			while (node_dof != node_dof_end) {
				*gnode_dof = *node_dof;
				++node_dof; ++gnode_dof;
			}
			++node_release;
		}
		++elem_release;
	}

	math::vector<NodeAllReleases> node_releases_vec(node_releases.size());
	
	auto node = node_releases.begin(), node_end = node_releases.end();
	auto node_vec = node_releases_vec.begin();
	while (node != node_end) {
		node->second.node_globalid = node->first;
		*node_vec = std::move(node->second);
		++node; ++node_vec;
	}

	std::sort(node_releases_vec.begin(),node_releases_vec.end());
	std::sort(elements_releases.begin(), elements_releases.end());
	
	return node_releases_vec;
}


void ModelTraits::parse_nodal_data(const Model& model, Assemble& assemble) {
    // matrix with nodes global dofs
    math::vector_t<size_t,2> nodes_dofs = math::zeros<size_t>(model.nodes_info.size(),BaseNode::DOFS);
    
    // iterators
    auto gnode = nodes_dofs.begin(), node_end = nodes_dofs.end();
	size_t gnodeID = 0;
    typename math::vector<size_t>::iterator gdof, gdof_end;

    auto node_info = model.nodes_info.begin();
    auto constraint_node = model.constraints_info.begin();
	typename decltype(NodeConstraint::dofs)::const_iterator constraint_dof;
    
	auto node_release = model.nodes_releases.begin()
		,node_release_end = model.nodes_releases.end();
	decltype(model.nodes_releases[0].dofs.begin()) released_dof,released_dof_end;
	// set nodes_dofs matrix by global dofs if current dof is not constrained
	// constrained dofs value = 0
    size_t global_id = 1;
    while (gnode != node_end) {
        gdof = gnode->begin();
		gdof_end = gnode->end();
        
		// if dof is constained, gdof=0
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

		// Released dof
		if ((node_release < node_release_end) && (gnodeID == node_release->node_globalid)) {
			// resize vector, iterator gnode not invalidate
			gnode->resize(BaseNode::DOFS + node_release->dofs.size());
			
			gdof = gnode->begin()+BaseNode::DOFS; // skip already desined dofs
			released_dof = node_release->dofs.begin();
			released_dof_end = node_release->dofs.end();
			// loop over released dof
			while (released_dof != released_dof_end) {
				
				if (*released_dof == NodeAllReleases::IsNotReleased) {
					*gdof = 0;
				} else {
					*gdof = global_id;
					++global_id;
				}
				++released_dof; ++gdof;
			}
			
			++node_release;
		}


		++gnode; ++node_info; ++gnodeID;
    }
	
	// Calculate elements properties
	auto elem = model.elements.begin(), elem_end = model.elements.end();
#ifdef BEAM_RELEASES
	auto elem_release = model.elements_releases.begin()
	    ,elem_release_end = model.elements_releases.end();
#endif

	for (;elem != elem_end; ++elem) {
		
#ifdef BEAM_RELEASES
		if ((elem_release < elem_release_end) && ((*elem)->nodes_dofs == elem_release->elem_id)) {
			(*elem)->calc_parameters(model.nodes_info,*elem_release);
			++elem_release;
			
			continue;
		}
#endif
		
		(*elem)->calc_parameters(model.nodes_info);
	}

	// Added dofs according to beam releases dofs
	// new dofs are added to back of the node already defined dofs
	// node_dofs = [g1,g2,g3,g4,g5,g6, r1,r2,r3,r4,r5,r6]
	// where gi - defined dofs
	// ri = 0 if dof is not released
	// ri = gdofID if dof is released
	// so element reliased dofs = gi + BaseNode::DOFS = gi + 6
	


	assemble.ndofs = global_id-1;
    assemble.nodes_dofs = std::move(nodes_dofs);
}


void ModelTraits::store_load_vector(const Model& model
								  , const Assemble& assemble
                                  , math::vector<double>& load)
{
	auto load_node     = model.loads_info.begin()
		,load_node_end = model.loads_info.end();
		
	decltype(model.loads_info[0].dofs.begin()) load_dof, load_dof_end;
	
	// decltype(GDofs.nodes_dofs.begin()) gnode;
	decltype(assemble.nodes_dofs[0].begin()) gdof;
	
	while (load_node < load_node_end) {
		load_dof = load_node->dofs.begin();
		load_dof_end = load_node->dofs.end();
		gdof = assemble.nodes_dofs[load_node->id].begin();
		while (load_dof < load_dof_end) {
			if (*gdof > 0) {
				load[*gdof-1] = *load_dof;
			}
			++load_dof; ++gdof;
		}
		++load_node;
	}
}
void ModelTraits::store_load_vector(Model& model
								  , const Assemble& assemble
                                  , math::vector<double>& load
								  , double time)
{
	auto load_node     = model.loads_var_info.begin()
		,load_node_end = model.loads_var_info.end();
		
	decltype(model.loads_var_info[0].dofs.begin()) load_dof, load_dof_end;
	
	// decltype(GDofs.nodes_dofs.begin()) gnode;
	decltype(assemble.nodes_dofs[0].begin()) gdof;
	
	while (load_node < load_node_end) {
		(*load_node)(time);
		load_dof = load_node->dofs.begin();
		load_dof_end = load_node->dofs.end();
		gdof = assemble.nodes_dofs[load_node->id].begin();
		while (load_dof < load_dof_end) {
			if (*gdof > 0) {
				load[*gdof-1] = *load_dof;
			}
			++load_dof; ++gdof;
		}
		++load_node;
	}
}


void ModelTraits::assemble_precomputing(const Model& model, Assemble& assemble) {
	/*
	Assemble global matrices by two loops
	1. define elements global dofs and global matrix column heights
	2. fill global matrix
	
	Store global matrix as its upper triangular part
	*/
	
	// Elements global dofs
	math::vector_t<size_t,2> ElemsGDofs(model.elements.size()); // might be not rectangular matrix!
	// Global matrix columns heights
	math::vector<size_t> colhs(assemble.ndofs);
	size_t min_dof;
	
	// 1-st Loop over all element
	auto elem = model.elements.begin(), elem_end = model.elements.end();
	auto elemdofs = ElemsGDofs.begin(); // iterator oncurrent element dofs
	decltype(ElemsGDofs[0].begin()) elemdof, elemdof_end; // iterator on current element dof
	decltype(model.elements[0]->nodes.begin()) elemnode, elemnode_end; // iterator on element node
	decltype(assemble.nodes_dofs[0].begin()) nodegdof, nodegdof_end, nodegdof_released; // iterator on element global dof
	
	auto elem_release = model.elements_releases.begin()
		,elem_release_end = model.elements_releases.end();
	decltype(model.elements_releases[0].nodes.begin()) node_release,node_release_end;
	size_t elemnode_id;
	size_t max_elem_dofs = 0;
	size_t elem_ndofs;
	// loop over all elements
	while (elem != elem_end) {
		*elemdofs = math::zeros<size_t>((*elem)->ndofs()); // init elemdofs array
		elemdof = elemdofs->begin();
		elemnode = (*elem)->nodes.begin();
		elemnode_end = (*elem)->nodes.end();
		// loop over element nodes
		while (elemnode != elemnode_end) {
			nodegdof = assemble.nodes_dofs[*elemnode].begin();
			nodegdof_end = nodegdof + (*elem)->ndofs_node();

			while (nodegdof != nodegdof_end) {
				*elemdof = *nodegdof;
				++elemdof;
				++nodegdof;
			}
			++elemnode;
		}

		// released dofs
		// std::cout << "elem nodes_dofs " << (*elem)->nodes_dofs << ", elem_release nodes_dofs " << elem_release->elem_id << std::endl;
		if ( (elem_release < elem_release_end) && ((*elem)->nodes_dofs == elem_release->elem_id)) {
			
			// loop over released element nodes
			node_release = elem_release->nodes.begin();
			node_release_end = elem_release->nodes.end();
			while (node_release != node_release_end) {
				// loop over released dofs
				elemdof = elemdofs->begin() + node_release->node_localid*BaseNode::DOFS; // start of dofs of current element node
				elemnode_id = (*elem)->nodes[node_release->node_localid];                // id of node with releases
				
				nodegdof_released = assemble.nodes_dofs[elemnode_id].begin() + (node_release->shift+1) * BaseNode::DOFS; // start of released gdofs
				nodegdof = assemble.nodes_dofs[elemnode_id].begin();													  // start of not released gdofs
				nodegdof_end = nodegdof + (*elem)->ndofs_node();
				// loop over all node dofs
				// std::cout << "\t" << GDofs.nodes_dofs[elemnode_id] << std::endl;
				while (nodegdof != nodegdof_end) {
					// if dof was released,  
					if (*nodegdof_released != 0) {
						*elemdof = *nodegdof_released;
					} 
					++nodegdof; ++nodegdof_released; ++elemdof;
				}
				
				++node_release;
			}
			
			++elem_release;
		}

		// calc current column heights
		min_dof = BaseElement::min_elem_dof(elemdofs->begin(),elemdofs->end());
		elemdof = elemdofs->begin();
		elemdof_end = elemdofs->end();
		while (elemdof != elemdof_end) {
			if (*elemdof+1 > min_dof) {
				colhs[*elemdof-1] = std::max(*elemdof+1 - min_dof,colhs[*elemdof-1]);
			}
			++elemdof;
		}

		// calc max element dofs
		elem_ndofs = (*elem)->ndofs();
		if (elem_ndofs > max_elem_dofs)
			max_elem_dofs = elem_ndofs;

		++elem; ++elemdofs;
	}


	math::vector<size_t> diags(assemble.ndofs+1);
	auto diag = diags.begin(), diag_end = diags.end()-1;
	auto colh_1 = colhs.begin();
	*diag = 0;
	++diag;
	while (diag != diag_end) {
		*diag = *(diag-1) + *colh_1;
		++diag;	++colh_1;
	}
	diags.last() = diags[assemble.ndofs-1] + *colh_1;


	// GDofs.stiffness = math::zeros<double>(sum(colhs));
	assemble.elems_dofs = std::move(ElemsGDofs);
	assemble.diags = std::move(diags);
	assemble.band_width = *std::max_element(colhs.begin(),colhs.end());
	assemble.colhs = std::move(colhs);
	
	
	assemble.max_elem_dofs = max_elem_dofs;
	assemble.elem_matrix = math::zeros<double>(max_elem_dofs,max_elem_dofs);
	assemble.elem_matrix2 = math::zeros<double>(max_elem_dofs,max_elem_dofs);
	assemble.elem_load = math::zeros<double>(max_elem_dofs);
	assemble.elem_displ = math::zeros<double>(max_elem_dofs);
	assemble.elem_vel = math::zeros<double>(max_elem_dofs);
	assemble.elem_accel = math::zeros<double>(max_elem_dofs);
	

}

/* Store element displacements vector from global displacements vector*/
void ModelTraits::store_element_state_vectors(typename math::vector<double>::iterator elem_displacement		/* Element displacements vector */
											, typename math::vector<size_t>::const_iterator elemgdof  		/* Element global dofs id*/
											, typename math::vector<size_t>::const_iterator elemgdof_end		/* Element global dofs boundary id*/
											, const    math::vector<double>& displacement)							/* Global displacements vector*/
{
	size_t dof;
	/* Loop over element dofs */
	for (;elemgdof < elemgdof_end
				;++elem_displacement	// next element dof value
				,++elemgdof)	// next element global dof id
	{
			dof = *elemgdof;
			// if dof is constrained - displacement == 0
			if (dof == Model::DOF_IS_CONSTRAINED) {
				*elem_displacement = 0.0;
				continue;					// go to next element dof
			}
			--dof;
			// if dof is not constrained
			*elem_displacement = displacement[dof];	// store element dof displacement value with `elemgdof` id
	}
}

void ModelTraits::store_element_state_vectors(typename math::vector<double>::iterator elem_displacement
											, typename math::vector<double>::iterator elem_velocity
								  			, typename math::vector<size_t>::const_iterator elemgdof
								  			, typename math::vector<size_t>::const_iterator elemgdof_end
								  			, const    math::vector<double>& displacement
											, const    math::vector<double>& velocity)
{
	size_t dof;
	/* Loop over element dofs */
	for (;elemgdof < elemgdof_end
				;++elem_displacement	// next element dof value
				,++elem_velocity
				,++elemgdof)	// next element global dof id
	{
			// if dof is constrained - displacement == 0
			dof = *elemgdof;
			if (dof == Model::DOF_IS_CONSTRAINED) {
				*elem_displacement = 0.0;
				*elem_velocity = 0.0;
				continue;					// go to next element dof
			}
			// if dof is not constrained
			--dof;
			*elem_displacement 	= displacement[dof];	// store element dof displacement value with `elemgdof` id
			*elem_velocity 		= velocity[dof];
	}
}

void ModelTraits::store_element_state_vectors(typename math::vector<double>::iterator elem_displacement
											, typename math::vector<double>::iterator elem_velocity
											, typename math::vector<double>::iterator elem_acceleration
								  			, typename math::vector<size_t>::const_iterator elemgdof
								  			, typename math::vector<size_t>::const_iterator elemgdof_end
								  			, const    math::vector<double>& displacement
											, const    math::vector<double>& velocity
											, const    math::vector<double>& acceleration)
{
	size_t dof;
	/* Loop over element dofs */
	for (;elemgdof < elemgdof_end
				;++elem_displacement	// next element dof value
				,++elem_velocity
				,++elem_acceleration
				,++elemgdof)	// next element global dof id
	{
			dof = *elemgdof;
			// if dof is constrained - displacement == 0
			if (dof == Model::DOF_IS_CONSTRAINED) {
				*elem_displacement = 0.0;
				*elem_velocity = 0.0;
				*elem_acceleration = 0.0;
				continue;					// go to next element dof
			}
			--dof;
			// if dof is not constrained
			*elem_displacement 	= displacement[dof];	// store element dof displacement value with `elemgdof` id
			*elem_velocity 		= velocity[dof];
			*elem_acceleration 	= acceleration[dof];
	}
}

} // namespace fem