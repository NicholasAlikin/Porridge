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
        ,elem_matrix2(other.elem_matrix2)
		,elem_load(other.elem_load)
        ,elem_load2(other.elem_load2)
		
        ,elem_displ(other.elem_displ)
        ,elem_displ2(other.elem_displ2)
        ,elem_vel(other.elem_vel)
        ,elem_accel(other.elem_accel)
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
        ,elem_matrix2(std::move(other.elem_matrix2))
		,elem_load(std::move(other.elem_load))
        ,elem_load2(std::move(other.elem_load2))
		
        ,elem_displ(std::move(other.elem_displ))
        ,elem_displ2(std::move(other.elem_displ2))
        ,elem_vel(std::move(other.elem_vel))
        ,elem_accel(std::move(other.elem_accel))
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


void ModelTraits::store_load_vector(Model& model
								  , const Assemble& assemble
                                  , math::vector<double>& load
								  , double time
                                  , const npath::DFT& dft)
{
	/* Iterator on loaded nodes */
    auto load_node     = model.loads_var_info.begin()
		,load_node_end = model.loads_var_info.end();
	
    /* Iterator on node load vector components */
	decltype(model.loads_var_info[0].dofs.begin()) load_dof, load_dof_end, load_hdof;
	
	/* Iterator on nodes dofs */
	decltype(assemble.nodes_dofs[0].begin()) gdof, gdof_end;
	
    /* Slices on load dof for current harmonic */
    decltype(load.begin()) gload;
    math::Slice sload(load.begin(), load.end(), assemble.ndofs);    // global
    math::Slice <decltype(model.loads_var_info[0].dofs.begin()),
                 decltype(model.loads_var_info[0].dofs.end())  > node_sload;  // element
    
    /* Loop over all loaded nodes */
	while (load_node < load_node_end) {
        /* Calculate load vector which corresponds (pseudo) `time` */
		(*load_node)(time);

		/* Node load slice end remains constant */
        // node_sload.update_to(load_node->dofs.end());
        /* Loop over node global dofs */
		for (gdof     = assemble.nodes_dofs[load_node->id].begin()
            ,gdof_end = assemble.nodes_dofs[load_node->id].end()
            ,load_dof = load_node->dofs.begin()
                    ;gdof < gdof_end
                            ;++gdof
                            ,++load_dof)
        {
            /* constraint dofs cannot be loaded */
            if (*gdof == Model::DOF_IS_CONSTRAINED) continue;
            /* update slices begins */
            /* Loop over harmonics */
            for (load_dof_end = load_node->dofs.end()
                ,load_hdof = load_dof
                ,gload = load.begin() + *gdof-1
                        ;load_hdof < load_dof_end
                                ;load_hdof += BaseNode::DOFS // step size
                                ,gload += assemble.ndofs) // step size
            {
                *gload = *load_hdof;
            }
            // sload.update_from(load.begin() + *gdof-1);
            // node_sload.new_slice(load_dof,load_node->dofs.end(),BaseNode::DOFS);
            // sload = node_sload; // not += because node can be loaded only by one load
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
	assemble.elem_matrix_size(max_elem_dofs);
    assemble.elem_load_size(max_elem_dofs);
    assemble.elem_state_size(max_elem_dofs);

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

void ModelTraits::store_element_state_vectors(typename math::vector<double>::iterator elem_displacement
											, typename math::vector<double>::iterator elem_velocity
											, typename math::vector<double>::iterator elem_acceleration
								  			, typename math::vector<size_t>::const_iterator elemgdof
								  			, typename math::vector<size_t>::const_iterator elemgdof_end
								  			, const    math::vector<double>& displacement
											, const    math::vector<double>& velocity
											, const    math::vector<double>& acceleration
                                            , const ::npath::DFT&            dft
                                            , size_t elem_ndofs)
{
	size_t dof;
    /* Slice to = Slice from for each element dof */
    // element state components
    size_t state_vec_size = elem_ndofs * dft.N;
    math::Slice elem_displ(elem_displacement, elem_displacement   +state_vec_size, elem_ndofs);
    math::Slice elem_vel(  elem_velocity,     elem_velocity       +state_vec_size, elem_ndofs);
    math::Slice elem_accel(elem_acceleration, elem_acceleration   +state_vec_size, elem_ndofs);
	// global state components
    math::Slice displ(displacement.begin(), displacement.end(), dft.ndof);
    math::Slice vel(  velocity.begin(),     velocity.end(),     dft.ndof);
    math::Slice accel(acceleration.begin(), acceleration.end(), dft.ndof);
	// std::cout << "elem_displ = " << elem_displ.size() << std::endl;
    // std::cout << "elem_vel = " << elem_vel.size() << std::endl;
    // std::cout << "elem_accel = " << elem_accel.size() << std::endl;
    // std::cout << "displ = " << displ.size() << std::endl;
    // std::cout << "vel = " << vel.size() << std::endl;
    // std::cout << "accel = " << accel.size() << std::endl;
    /* Loop over element dofs */
	for (;elemgdof < elemgdof_end
				;++elem_displacement	// next element dof value
				,++elem_velocity
				,++elem_acceleration
				,++elemgdof)	// next element global dof id
	{
			dof = *elemgdof;

            elem_displ.update_from(elem_displacement);
            elem_vel.update_from(elem_velocity);
            elem_accel.update_from(elem_acceleration);
			
            // if dof is constrained - displacement == 0
			if (dof == Model::DOF_IS_CONSTRAINED) {
				elem_displ = 0.0;
                elem_vel = 0.0;
                elem_accel = 0.0;
				continue;					// go to next element dof
			}
			--dof;
            displ.update_from(displacement.begin() +dof);
            vel.update_from(velocity.begin()       +dof);
            accel.update_from(acceleration.begin() +dof);
			// if dof is not constrained
			elem_displ 	= displ;	// store element dof displacement value with `elemgdof` id
			elem_vel    = vel;
			elem_accel 	= accel;
	}
}


void ModelTraits::store_element_state_vectors(typename math::vector<double>::iterator elem_displacement
											, typename math::vector<double>::iterator elem_velocity
											, typename math::vector<double>::iterator elem_acceleration
                                            , typename math::vector<double>::iterator elem_displ_freq
								  			, typename math::vector<size_t>::const_iterator elemgdof
								  			, typename math::vector<size_t>::const_iterator elemgdof_end
								  			, const    math::vector<double>& displacement
											, const    math::vector<double>& velocity
											, const    math::vector<double>& acceleration
                                            , const    math::vector_const_slice<double>& displ_freq
                                            , const ::npath::DFT& dft
                                            , size_t elem_ndofs)
{
	size_t dof;
    /* Slice to = Slice from for each element dof */
    // element state components
    size_t state_vec_size       = elem_ndofs * dft.time_basic_size();
    size_t state_vec_freq_size  = elem_ndofs * dft.frequency_basic_size();

    math::Slice elem_displ(elem_displacement,
                           elem_displacement + state_vec_size, elem_ndofs);
    
    math::Slice elem_vel(  elem_velocity,     
                           elem_velocity + state_vec_size, elem_ndofs);
    
    math::Slice elem_accel(elem_acceleration, 
                           elem_acceleration + state_vec_size, elem_ndofs);
    
    math::Slice elem_displfreq(elem_displ_freq,
                               elem_displ_freq +state_vec_freq_size, elem_ndofs);

	// global state components
    math::Slice displ(displacement.begin(), displacement.end(), dft.ndof);
    math::Slice vel(  velocity.begin(),     velocity.end(),     dft.ndof);
    math::Slice accel(acceleration.begin(), acceleration.end(), dft.ndof);
	math::Slice displfreq(displ_freq.begin(), displ_freq.end(), dft.ndof);
	
    /* Loop over element dofs */
	for (;elemgdof < elemgdof_end
				;++elem_displacement	// next element dof value
				,++elem_velocity
				,++elem_acceleration
                ,++elem_displ_freq
				,++elemgdof)	// next element global dof id
	{
			dof = *elemgdof;

            elem_displ.update_from(elem_displacement);
            elem_vel.update_from(elem_velocity);
            elem_accel.update_from(elem_acceleration);
            elem_displfreq.update_from(elem_displ_freq);
			
            // if dof is constrained - displacement == 0
			if (dof == Model::DOF_IS_CONSTRAINED) {
				elem_displ = 0.0;
                elem_vel = 0.0;
                elem_accel = 0.0;
                elem_displfreq = 0.0;
				continue;					// go to next element dof
			}
			--dof;
            displ.update_from(displacement.begin()  +dof);
            vel.update_from(velocity.begin()        +dof);
            accel.update_from(acceleration.begin()  +dof);
            displfreq.update_from(displ_freq.begin()+dof);
			// if dof is not constrained
			elem_displ 	= displ;	// store element dof displacement value with `elemgdof` id
			elem_vel    = vel;
			elem_accel 	= accel;
            elem_displfreq = displfreq;
	}
}



/* Assemble global tangent matrix and load vector in frequency domain. */
void ModelTraits::assemble(const Model& model, const Assemble& assemble
                                ,       math::vector_t<double,2>&   matrix              /* tangent matrix in frequency domain*/
                                ,       math::vector_t<double,1>&   load                /* load vector in frequency domain */
                                , const math::vector_t<double,1>&   q                   /* displacements in time domain */
			                    , const math::vector_t<double,1>&   dqdt                /* velocities in time domain */
								, const math::vector_t<double,1>&   d2qdt2              /* accelerations in time domain */
                                , const ::npath::DFT&               dft                 /* DFT object */
                                , double                            freq                /* current frequency */
                                ,       math::vector_t<double,3>&   buffer_dft_matrix   /* buffer to calculate DFT of matrix */
                                ,       math::vector_t<double,1>&   buffer_dft_vector   /* buffer to calculate DFT of vector */
                                , void (BaseElement::* element_matrix_load)(            /* Calculate element local matrix and load vector: */
                                                       math::vector_t<double,2>&        /*      element local matrix */
												,      math::vector_t<double,1>&        /*      element local load */ 
												,const math::vector_t<double,1>&        /*      element property */ 
												,const math::vector_t<double,1>&        /*      element material */
												,const math::vector_t<double,1>&        /*      element q */
												,const math::vector_t<double,1>&        /*      element dqdt */
												,const math::vector_t<double,1>&        /*      element d2qdt2 */  
                                                ,const ::npath::DFT&                    /*      dft */ 
                                                ,double                                 /*      freq */ 
                                                ,      math::vector_t<double,3>&        /*      buffer_dft_matrix */
                                                ,      math::vector_t<double,1>&        /*      buffer_dft_vector */) const
                                )
{
	/*Iterators*/	
	decltype(assemble.elem_load.begin())  		load_dof;								/* Element load vector */
 
	decltype(assemble.elem_matrix.begin()) 		matloc_row, matloc_hrow;				/* Element matrix row values: inital (zeros harmonic) and current harmonic */
	decltype(assemble.elem_matrix[0].begin()) 	matloc_col;				                /* Element matrix col values */
    decltype(matrix.begin())                    matrix_hrow;                            /* Global  matrix row current harmonic values */  
    decltype(matrix[0].begin())                 matrix_col;                             /* Global  matrix row values */

	auto 										elem = model.elements.begin()			/* Element*/
											  , elem_end = model.elements.end();	
	
	auto 										elemgdofs = assemble.elems_dofs.begin();/* Element dofs */
	decltype(assemble.elems_dofs[0].begin()) 	elemgdof_row, elemgdof_row_end 			/* Element dof for loop over element matrix rows */
											  , elemgdof_col, elemgdof_col_end;  		/*             for loop over element matrix cols */
	
    /*Slices. In time domain assemble needs loop over dofs,
    in frequency domain - loop over harmonics for each dofs is needed also.
    `s` - means Slice 
    This slices containts values in frequency domain of certain dof and all harmonics*/
    // const slices step, equal to global ndofs
    size_t hndofs = dft.frequency_size();
    math::Slice sload(  load.begin(), load.begin()+hndofs, assemble.ndofs);
    math::Slice smatrix(matrix[0].begin(), matrix[0].begin()+hndofs,
                        assemble.ndofs);
    
    // step equal to element ndofs
    math::Slice<decltype(assemble.elem_load.begin()),
                decltype(assemble.elem_load.end())>         elem_sload;
    math::Slice<decltype(assemble.elem_matrix[0].begin()),
                decltype(assemble.elem_matrix[0].end())>    elem_smatrix;
    /* Other variables */
    size_t elem_ndofs;

	/* Loop over all elements */
	while (elem != elem_end) {
        elem_ndofs = (*elem)->ndofs();
        /* store element state vectors */
		ModelTraits::store_element_state_vectors(assemble.elem_displ.begin()	// where store to
												,assemble.elem_vel.begin()
												,assemble.elem_accel.begin()
												,elemgdofs->begin()				// with dofs id
												,elemgdofs->end()
												,q								// store from there
												,dqdt
												,d2qdt2
                                                ,dft
                                                ,elem_ndofs);
        /* calculate element matrix and load vector in frequency domain */
		((*elem)->*element_matrix_load)( assemble.elem_matrix 					// where to store matrix
										,assemble.elem_load						// where to store load vector
										,model.properties[(*elem)->propID]		// element property
										,model.materials[ (*elem)->matlID]		// element material
										,assemble.elem_displ					// element displacement vector
										,assemble.elem_vel
										,assemble.elem_accel
                                        ,dft
                                        ,freq
                                        ,buffer_dft_matrix
                                        ,buffer_dft_vector);
              
		/* store element matrix and load vector to global matrix and load vector in frequency domain 
         like it in time domain */
        // update element matrix and load slices
        // - because in inner loops this slices remain const step,
        // equal to element ndofs
        elem_sload.new_slice(assemble.elem_load.begin()
                            ,assemble.elem_load.begin()
                                 +elem_ndofs*dft.frequency_basic_size()
                            ,elem_ndofs);
        elem_smatrix.new_slice(assemble.elem_matrix[0].begin(),
                               assemble.elem_matrix[0].begin()
                                   +elem_ndofs*dft.frequency_basic_size(),
                               elem_ndofs);
		/* loop over inital matrix rows (like in time domain) */
		for (matloc_row         = assemble.elem_matrix.begin()    	// iterator on elem matrix row
		    ,elemgdof_row       = elemgdofs->begin()   				// iterator on global dofs, corresponding row global dofs of the elem
			,elemgdof_row_end   = elemgdofs->end() 					// same
			,load_dof           = assemble.elem_load.begin()        // iterator on elem load vector
		   				;elemgdof_row < elemgdof_row_end 			// Loop over all matrix rows
									;++matloc_row
									,++elemgdof_row
									,++load_dof)
		{

			// if dof is constrained - go to next matrix row and vector component
			if (*elemgdof_row == Model::DOF_IS_CONSTRAINED) continue;     					// assamble only matrix rows, which correspond not constrained global dofs
            
            /* Loop over harmonics rows */
            for (size_t h = 0; h < dft.frequency_basic_size(); ++h) {
                matloc_hrow = matloc_row + h*elem_ndofs;
                matrix_hrow = matrix.begin() + *elemgdof_row-1 + h*assemble.ndofs;

                // update matrix Slices end bounds
                matrix_col = matrix_hrow->begin(); //matrix[*elemgdof_row-1].begin();
                smatrix.update_to(matrix_col + hndofs);
                elem_smatrix.update_to( matloc_hrow->end());
                // loop over column of the current matrix row
                // loop over all element - assume that matrix is not symmetric
                for (matloc_col = matloc_hrow->begin()  					// iterator on column element
                    ,elemgdof_col = elemgdofs->begin()					// iterator on global dofs, corresponding col global dofs of the elem
                    ,elemgdof_col_end = elemgdofs->end() 				// same
                                ;elemgdof_col < elemgdof_col_end 		// loop over row elements from diagonal to the end
                                            ;++matloc_col
                                            ,++elemgdof_col)
                {
                    // if dof is constrained - go to next component
                    if (*elemgdof_col == Model::DOF_IS_CONSTRAINED) continue;  	// assamble only matrix columns, which correspond not constrained global dofs
                    
                    // update matrix Slices
                    smatrix.update_from(matrix_col + *elemgdof_col-1);
                    elem_smatrix.update_from(matloc_col);
                    // assemble global matrix for each harmonic
                    smatrix += elem_smatrix;
                    // std::cout << "row " << *elemgdof_row << ", col " << *elemgdof_col << "elem_smatrix = " << elem_smatrix.size() << ", smatrix = " << smatrix.size() << std::endl;
                
                } // loop over column of the current matrix row

            } // Loop over harmonics rows

			/* store global load vector */
            // update load Slices
            sload.update_from(load.begin() + *elemgdof_row-1); // *elemgdof_row-1 == current dof number
            elem_sload.update_from(load_dof);
            // assemble global load for each harmonic
            sload += elem_sload;
        
        } // loop over inital matrix rows
		
		
		++elem; ++elemgdofs;

	} // Loop over all elements
	
}


/* Assemble extendent Jacobi matrix */
void ModelTraits::assemble(const Model& model, const Assemble& assemble
                            ,       math::vector_t<double,2>&   matrix      /* global system Jacobi matrix */
                            ,       math::vector_t<double,1>&   load        /* global internal system load vector */
                            , const math::vector_t<double,1>&   u           /* time domain displacement */
                            , const math::vector_t<double,1>&   dudt        /* time domain velocity */
                            , const math::vector_t<double,1>&   d2udt2      /* time domain acceleration */
                            , const math::vector_const_slice<double>& q           /* frequency domain displacement */
                            , const ::npath::DFT&               dft         /* DFT transformer */
                            , double                            freq        /* current frequency */
                            ,       math::vector_t<double,3>&   buffer_dft_matrix
                            ,       math::vector_t<double,1>&   buffer_dft_vector
                            , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                    math::vector_t<double,2>&   /* element local matrix */
                                            ,      math::vector_t<double,1>&    /* element local load */
                                            ,      math::vector_t<double,1>&    /* element extendent matrix column */ 
                                            ,const math::vector_t<double,1>&    /* element property */ 
                                            ,const math::vector_t<double,1>&    /* element material */
                                            ,const math::vector_t<double,1>&    /* element u */
                                            ,const math::vector_t<double,1>&    /* element dudt */
                                            ,const math::vector_t<double,1>&    /* element d2udt2 */
                                            ,const math::vector_const_slice<double>& /* element q */  
                                            ,const ::npath::DFT&                /* dft */ 
                                            ,double                             /* freq */ 
                                            ,      math::vector_t<double,3>&    /* buffer_dft_matrix */
                                            ,      math::vector_t<double,1>&    /* buffer_dft_vector */) const
                            )
{
	/*Iterators*/	
	decltype(assemble.elem_load.begin())  		load_dof;								/* Element load vector */
    decltype(assemble.elem_load2.begin())  		load2_dof,load2_hdof;					/* Element load2 vector current dof and corresponding harmonic */
    
 
	decltype(assemble.elem_matrix.begin()) 		matloc_row, matloc_hrow;				/* Element matrix row values: inital (zeros harmonic) and current harmonic */
	decltype(assemble.elem_matrix[0].begin()) 	matloc_col;				                /* Element matrix col values */
    decltype(matrix.begin())                    matrix_hrow;                            /* Global  matrix row current harmonic values */  
    decltype(matrix[0].begin())                 matrix_col;                             /* Global  matrix row values */

	auto 										elem = model.elements.begin()			/* Element*/
											  , elem_end = model.elements.end();	
	
	auto 										elemgdofs = assemble.elems_dofs.begin();/* Element dofs */
	decltype(assemble.elems_dofs[0].begin()) 	elemgdof_row, elemgdof_row_end 			/* Element dof for loop over element matrix rows */
											  , elemgdof_col, elemgdof_col_end;  		/*             for loop over element matrix cols */
	
    /*Slices. In time domain assemble needs loop over dofs,
    in frequency domain - loop over harmonics for each dofs is needed also.
    `s` - means Slice 
    This slices containts values in frequency domain of certain dof and all harmonics*/
    // const slices step, equal to global ndofs
    size_t hndofs = dft.frequency_size();
    math::Slice sload(  load.begin(), load.begin()+hndofs, assemble.ndofs);
    math::Slice smatrix(matrix[0].begin(), matrix[0].begin()+hndofs,
                        assemble.ndofs);
    
    // step equal to element ndofs
    math::Slice<decltype(assemble.elem_load.begin()),
                decltype(assemble.elem_load.end())>         elem_sload;
    math::Slice<decltype(assemble.elem_matrix[0].begin()),
                decltype(assemble.elem_matrix[0].end())>    elem_smatrix;
    /* Other variables */
    size_t elem_ndofs, elem_hndofs;

	/* Loop over all elements */
	while (elem != elem_end) {
        elem_ndofs = (*elem)->ndofs();
        elem_hndofs = elem_ndofs * dft.frequency_basic_size();
		/* store element state vectors */
        ModelTraits::store_element_state_vectors(assemble.elem_displ.begin()	// where store to
												,assemble.elem_vel.begin()
												,assemble.elem_accel.begin()
                                                ,assemble.elem_displ2.begin()
												,elemgdofs->begin()				// with dofs id
												,elemgdofs->end()
												,u								// store from there
												,dudt
												,d2udt2
                                                ,q
                                                ,dft
                                                ,elem_ndofs);
        math::Slice elem_hdispl(assemble.elem_displ2.begin(),
                                assemble.elem_displ2.begin()+elem_hndofs);
        /* calculate element matrix and load vector in frequency domain */
		((*elem)->*element_matrix_load)( assemble.elem_matrix 					// where to store matrix
										,assemble.elem_load						// where to store load vector
                                        ,assemble.elem_load2
										,model.properties[(*elem)->propID]		// element property
										,model.materials[ (*elem)->matlID]		// element material
										,assemble.elem_displ					// element displacement vector
										,assemble.elem_vel
										,assemble.elem_accel
                                        ,elem_hdispl
                                        ,dft
                                        ,freq
                                        ,buffer_dft_matrix
                                        ,buffer_dft_vector);
#if 0
        if ((*elem)->nodes_dofs == 0) {
            // std::cout << "1st elem: J = \n" << assemble.elem_matrix << '\n';
            // << "\ndr/dw = " << assemble.elem_load2 << '\n';

            /* nnumerical tangent matrix calculation */
            auto num_matrix = math::zeros<double>(assemble.elem_matrix);
            auto num_load = math::zeros<double>(assemble.elem_load);
            auto num_load2 = math::zeros<double>(assemble.elem_load);

            dft.time_domain(assemble.elem_displ2,freq,assemble.elem_displ,assemble.elem_vel,assemble.elem_accel,elem_ndofs);
            (*elem)->frequency_Load( num_load						// where to store load vector
                                    ,model.properties[(*elem)->propID]		// element property
                                    ,model.materials[ (*elem)->matlID]		// element material	
                                    ,assemble.elem_displ					// element displacement vector
                                    ,assemble.elem_vel
                                    ,assemble.elem_accel
                                    ,dft
                                    ,freq
                                    ,buffer_dft_vector);
            double dx = 1e-5;
            size_t dof;
            for (dof = 0; dof < assemble.elem_displ2.size(); ++dof) {
                assemble.elem_displ2[dof] += dx;
                dft.time_domain(assemble.elem_displ2,freq,assemble.elem_displ,assemble.elem_vel,assemble.elem_accel,elem_ndofs);
                math::fill(num_load2,0.0);
                (*elem)->frequency_Load( num_load2						// where to store load vector
                                    ,model.properties[(*elem)->propID]		// element property
                                    ,model.materials[ (*elem)->matlID]		// element material	
                                    ,assemble.elem_displ					// element displacement vector
                                    ,assemble.elem_vel
                                    ,assemble.elem_accel
                                    ,dft
                                    ,freq
                                    ,buffer_dft_vector);
                for (size_t i = 0; i < assemble.elem_displ2.size(); ++i) {
                    num_matrix[i][dof] = (num_load2[i]-num_load[i])/dx;
                }
                assemble.elem_displ2[dof] -= dx;
            }
            // dr/dw
            dft.time_domain(assemble.elem_displ2,freq+dx,assemble.elem_displ,assemble.elem_vel,assemble.elem_accel,elem_ndofs);
            math::fill(num_load2,0.0);
            (*elem)->frequency_Load( num_load2						// where to store load vector
                                    ,model.properties[(*elem)->propID]		// element property
                                    ,model.materials[ (*elem)->matlID]		// element material	
                                    ,assemble.elem_displ					// element displacement vector
                                    ,assemble.elem_vel
                                    ,assemble.elem_accel
                                    ,dft
                                    ,freq+dx
                                    ,buffer_dft_vector);
            auto num_drdw = (num_load2-num_load)/dx;

            // std::cout << "num J = \n" << num_matrix << '\n'
            
            std::cout << "# |dr/dw - num dr/dw| = " << math::norm(assemble.elem_load2 - num_drdw)/math::norm(num_drdw)
            << ", |dr/dw| = " << math::norm(assemble.elem_load2)
            << ", |num dr/dw| = " << math::norm(num_drdw) << std::endl;
            // << "\nnum dr/dw = " << num_drdw << std::endl;
            assemble.elem_load2 = num_drdw;
            assemble.elem_matrix = num_matrix;
        }
#endif
		/* store element matrix and load vector to global matrix and load vector in frequency domain 
         like it in time domain */
        // update element matrix and load slices
        // - because in inner loops this slices remain const step,
        // equal to element ndofs
        elem_sload.new_slice(assemble.elem_load.begin()
                            ,assemble.elem_load.begin()
                                 +elem_ndofs*dft.frequency_basic_size()
                            ,elem_ndofs);
        elem_smatrix.new_slice(assemble.elem_matrix[0].begin(),
                               assemble.elem_matrix[0].begin()
                                   +elem_ndofs*dft.frequency_basic_size(),
                               elem_ndofs);
		/* loop over inital matrix rows (like in time domain) */
		for (matloc_row         = assemble.elem_matrix.begin()    	// iterator on elem matrix row
		    ,elemgdof_row       = elemgdofs->begin()   				// iterator on global dofs, corresponding row global dofs of the elem
			,elemgdof_row_end   = elemgdofs->end() 					// same
			,load_dof           = assemble.elem_load.begin()        // iterator on elem load vector
            ,load2_dof          = assemble.elem_load2.begin()       // iterator on elem load2 vector (here load2 is extendent matrix column)
		   				;elemgdof_row < elemgdof_row_end 			// Loop over all matrix rows
									;++matloc_row
									,++elemgdof_row
									,++load_dof
                                    ,++load2_dof)
		{

			// if dof is constrained - go to next matrix row and vector component
			if (*elemgdof_row == Model::DOF_IS_CONSTRAINED) continue;     					// assamble only matrix rows, which correspond not constrained global dofs
            
            /* Loop over harmonics rows */
            for (size_t h = 0; h < dft.frequency_basic_size(); ++h) {
                matloc_hrow = matloc_row + h*elem_ndofs;
                matrix_hrow = matrix.begin() + *elemgdof_row-1 + h*assemble.ndofs;
                
                load2_hdof = load2_dof + h*elem_ndofs;

                // update matrix Slices end bounds
                matrix_col = matrix_hrow->begin(); //matrix[*elemgdof_row-1].begin();
                smatrix.update_to(matrix_col + hndofs);
                elem_smatrix.update_to(matloc_hrow->end());
                // loop over column of the current matrix row
                // loop over all element - assume that matrix is not symmetric
                for (matloc_col = matloc_hrow->begin()  					// iterator on column element
                    ,elemgdof_col = elemgdofs->begin()					// iterator on global dofs, corresponding col global dofs of the elem
                    ,elemgdof_col_end = elemgdofs->end() 				// same
                                ;elemgdof_col < elemgdof_col_end 		// loop over row elements from diagonal to the end
                                            ;++matloc_col
                                            ,++elemgdof_col)
                {
                    // if dof is constrained - go to next component
                    if (*elemgdof_col == Model::DOF_IS_CONSTRAINED) continue;  	// assamble only matrix columns, which correspond not constrained global dofs
                    
                    // update matrix Slices
                    smatrix.update_from(matrix_col + *elemgdof_col-1);
                    elem_smatrix.update_from(matloc_col);
                    // assemble global matrix for each harmonic
                    smatrix += elem_smatrix;

                
                } // loop over column of the current matrix row
                
                /* store extendent column */
                matrix_hrow->last() += *load2_hdof;

            } // Loop over harmonics rows

			/* store global load vector */
            // update load Slices
            sload.update_from(load.begin() + *elemgdof_row-1); // *elemgdof_row-1 == current dof number
            elem_sload.update_from(load_dof);
            // assemble global load for each harmonic
            sload += elem_sload;
        
        } // loop over inital matrix rows
		// if ((*elem)->nodes_dofs == 0) {
        //     std::cout << "\n...1st elem: Jglobal = \n" << matrix << '\n';
        // }
		
		++elem; ++elemgdofs;

	} // Loop over all elements
	
}


/* Assemble load only */
void ModelTraits::assemble(const Model& model, const Assemble& assemble
                            ,       math::vector_t<double,1>&   load        /* global internal system load vector */
                            , const math::vector_t<double,1>&   u           /* time domain displacement */
                            , const math::vector_t<double,1>&   dudt        /* time domain velocity */
                            , const math::vector_t<double,1>&   d2udt2      /* time domain acceleration */
                            , const ::npath::DFT&               dft         /* DFT transformer */
                            , double                            freq        /* current frequency */
                            ,       math::vector_t<double,1>&   buffer_dft_vector
                            , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                   math::vector_t<double,1>&    /* element local load */
                                            ,const math::vector_t<double,1>&    /* element property */ 
                                            ,const math::vector_t<double,1>&    /* element material */
                                            ,const math::vector_t<double,1>&    /* element u */
                                            ,const math::vector_t<double,1>&    /* element dudt */
                                            ,const math::vector_t<double,1>&    /* element d2udt2 */
                                            ,const ::npath::DFT&                /* dft */ 
                                            ,double                             /* freq */ 
                                            ,      math::vector_t<double,1>&    /* buffer_dft_vector */) const
                            )
{
	/*Iterators*/	
	decltype(assemble.elem_load.begin())  		load_dof;								/* Element load vector */
    decltype(assemble.elem_load2.begin())  		load2_dof,load2_hdof;					/* Element load2 vector current dof and corresponding harmonic */
    
 
	auto 										elem = model.elements.begin()			/* Element*/
											  , elem_end = model.elements.end();	
	
	auto 										elemgdofs = assemble.elems_dofs.begin();/* Element dofs */
	decltype(assemble.elems_dofs[0].begin()) 	elemgdof_row, elemgdof_row_end 			/* Element dof for loop over element matrix rows */
											  , elemgdof_col, elemgdof_col_end;  		/*             for loop over element matrix cols */
	
    /*Slices. In time domain assemble needs loop over dofs,
    in frequency domain - loop over harmonics for each dofs is needed also.
    `s` - means Slice 
    This slices containts values in frequency domain of certain dof and all harmonics*/
    // const slices step, equal to global ndofs
    size_t hndofs = dft.frequency_size();
    math::Slice sload(  load.begin(), load.begin()+hndofs, assemble.ndofs);
    
    // step equal to element ndofs
    math::Slice<decltype(assemble.elem_load.begin()),
                decltype(assemble.elem_load.end())>         elem_sload;
    
    /* Other variables */
    size_t elem_ndofs, elem_hndofs;

	/* Loop over all elements */
	while (elem != elem_end) {
        elem_ndofs = (*elem)->ndofs();
        elem_hndofs = elem_ndofs * dft.frequency_basic_size();
		/* store element state vectors */
        ModelTraits::store_element_state_vectors(assemble.elem_displ.begin()	// where store to
												,assemble.elem_vel.begin()
												,assemble.elem_accel.begin()
 												,elemgdofs->begin()				// with dofs id
												,elemgdofs->end()
												,u								// store from there
												,dudt
												,d2udt2
                                                ,dft
                                                ,elem_ndofs);
        math::Slice elem_hdispl(assemble.elem_displ2.begin(),
                                assemble.elem_displ2.begin()+elem_hndofs);
        /* calculate element matrix and load vector in frequency domain */
		((*elem)->*element_matrix_load)( assemble.elem_load
										,model.properties[(*elem)->propID]		// element property
										,model.materials[ (*elem)->matlID]		// element material
										,assemble.elem_displ					// element displacement vector
										,assemble.elem_vel
										,assemble.elem_accel
                                        ,dft
                                        ,freq
                                        ,buffer_dft_vector);
              
		/* store element matrix and load vector to global matrix and load vector in frequency domain 
         like it in time domain */
        // update element matrix and load slices
        // - because in inner loops this slices remain const step,
        // equal to element ndofs
        elem_sload.new_slice(assemble.elem_load.begin()
                            ,assemble.elem_load.begin()
                                 +elem_ndofs*dft.frequency_basic_size()
                            ,elem_ndofs);
        /* loop over inital matrix rows (like in time domain) */
		for (elemgdof_row       = elemgdofs->begin()   				// iterator on global dofs, corresponding row global dofs of the elem
			,elemgdof_row_end   = elemgdofs->end() 					// same
			,load_dof           = assemble.elem_load.begin()        // iterator on elem load vector
            			;elemgdof_row < elemgdof_row_end 			// Loop over all matrix rows
									;++elemgdof_row
									,++load_dof)
		{

			// if dof is constrained - go to next matrix row and vector component
			if (*elemgdof_row == Model::DOF_IS_CONSTRAINED) continue;     					// assamble only matrix rows, which correspond not constrained global dofs
            
			/* store global load vector */
            // update load Slices
            sload.update_from(load.begin() + *elemgdof_row-1); // *elemgdof_row-1 == current dof number
            elem_sload.update_from(load_dof);
            // assemble global load for each harmonic
            sload += elem_sload;
        
        } // loop over inital matrix rows
		
		
		++elem; ++elemgdofs;

	} // Loop over all elements
	
}


/* Assemble Jacobi as Eigen::SparseMatrix */
void ModelTraits::assemble(const Model& model, const Assemble& assemble
                                ,       Eigen::SparseMatrix<double>&   matrix              /* tangent matrix in frequency domain*/
                                ,       math::vector_t<double,1>&   load                /* load vector in frequency domain */
                                , const math::vector_t<double,1>&   q                   /* displacements in time domain */
			                    , const math::vector_t<double,1>&   dqdt                /* velocities in time domain */
								, const math::vector_t<double,1>&   d2qdt2              /* accelerations in time domain */
                                , const ::npath::DFT&               dft                 /* DFT object */
                                , double                            freq                /* current frequency */
                                ,       math::vector_t<double,3>&   buffer_dft_matrix   /* buffer to calculate DFT of matrix */
                                ,       math::vector_t<double,1>&   buffer_dft_vector   /* buffer to calculate DFT of vector */
                                , void (BaseElement::* element_matrix_load)(            /* Calculate element local matrix and load vector: */
                                                       math::vector_t<double,2>&        /*      element local matrix */
												,      math::vector_t<double,1>&        /*      element local load */ 
												,const math::vector_t<double,1>&        /*      element property */ 
												,const math::vector_t<double,1>&        /*      element material */
												,const math::vector_t<double,1>&        /*      element q */
												,const math::vector_t<double,1>&        /*      element dqdt */
												,const math::vector_t<double,1>&        /*      element d2qdt2 */  
                                                ,const ::npath::DFT&                    /*      dft */ 
                                                ,double                                 /*      freq */ 
                                                ,      math::vector_t<double,3>&        /*      buffer_dft_matrix */
                                                ,      math::vector_t<double,1>&        /*      buffer_dft_vector */) const
                                )
{
	/*Iterators*/	
	decltype(assemble.elem_load.begin())  		load_dof;								/* Element load vector */
 
	decltype(assemble.elem_matrix.begin()) 		matloc_row, matloc_hrow;				/* Element matrix row values: inital (zeros harmonic) and current harmonic */
	decltype(assemble.elem_matrix[0].begin()) 	matloc_col;				                /* Element matrix col values */
    
	auto 										elem = model.elements.begin()			/* Element*/
											  , elem_end = model.elements.end();	
	
	auto 										elemgdofs = assemble.elems_dofs.begin();/* Element dofs */
	decltype(assemble.elems_dofs[0].begin()) 	elemgdof_row, elemgdof_row_end 			/* Element dof for loop over element matrix rows */
											  , elemgdof_col, elemgdof_col_end;  		/*             for loop over element matrix cols */
	
    /*Slices. In time domain assemble needs loop over dofs,
    in frequency domain - loop over harmonics for each dofs is needed also.
    `s` - means Slice 
    This slices containts values in frequency domain of certain dof and all harmonics*/
    // const slices step, equal to global ndofs
    size_t hndofs = dft.frequency_size();
    math::Slice sload(  load.begin(), load.begin()+hndofs, assemble.ndofs);
    
    // step equal to element ndofs
    math::vector_slice<double> elem_sload;
    math::vector_slice<double> elem_smatrix;
    /* Other variables */
    size_t elem_ndofs, matrix_row;

	/* Loop over all elements */
	while (elem != elem_end) {
        elem_ndofs = (*elem)->ndofs();
        /* store element state vectors */
		ModelTraits::store_element_state_vectors(assemble.elem_displ.begin()	// where store to
												,assemble.elem_vel.begin()
												,assemble.elem_accel.begin()
												,elemgdofs->begin()				// with dofs id
												,elemgdofs->end()
												,q								// store from there
												,dqdt
												,d2qdt2
                                                ,dft
                                                ,elem_ndofs);
        /* calculate element matrix and load vector in frequency domain */
		((*elem)->*element_matrix_load)( assemble.elem_matrix 					// where to store matrix
										,assemble.elem_load						// where to store load vector
										,model.properties[(*elem)->propID]		// element property
										,model.materials[ (*elem)->matlID]		// element material
										,assemble.elem_displ					// element displacement vector
										,assemble.elem_vel
										,assemble.elem_accel
                                        ,dft
                                        ,freq
                                        ,buffer_dft_matrix
                                        ,buffer_dft_vector);
              
		/* store element matrix and load vector to global matrix and load vector in frequency domain 
         like it in time domain */
        // update element matrix and load slices
        // - because in inner loops this slices remain const step,
        // equal to element ndofs
        elem_sload.new_slice(assemble.elem_load.begin()
                            ,assemble.elem_load.begin()
                                 +elem_ndofs*dft.frequency_basic_size()
                            ,elem_ndofs);
        elem_smatrix.new_slice(assemble.elem_matrix[0].begin(),
                               assemble.elem_matrix[0].begin()
                                   +elem_ndofs*dft.frequency_basic_size(),
                               elem_ndofs);
		/* loop over inital matrix rows (like in time domain) */
		for (matloc_row         = assemble.elem_matrix.begin()    	// iterator on elem matrix row
		    ,elemgdof_row       = elemgdofs->begin()   				// iterator on global dofs, corresponding row global dofs of the elem
			,elemgdof_row_end   = elemgdofs->end() 					// same
			,load_dof           = assemble.elem_load.begin()        // iterator on elem load vector
		   				;elemgdof_row < elemgdof_row_end 			// Loop over all matrix rows
									;++matloc_row
									,++elemgdof_row
									,++load_dof)
		{

			// if dof is constrained - go to next matrix row and vector component
			if (*elemgdof_row == Model::DOF_IS_CONSTRAINED) continue;     					// assamble only matrix rows, which correspond not constrained global dofs
            
            /* Loop over harmonics rows */
            for (size_t h = 0; h < dft.frequency_basic_size(); ++h) {
                matloc_hrow = matloc_row + h*elem_ndofs;
                matrix_row = *elemgdof_row-1 + h*assemble.ndofs;
                // update matrix Slices end bounds
                elem_smatrix.update_to( matloc_hrow->end());
                // loop over column of the current matrix row
                // loop over all element - assume that matrix is not symmetric
                for (matloc_col = matloc_hrow->begin()  			    // iterator on column element
                    ,elemgdof_col = elemgdofs->begin()					// iterator on global dofs, corresponding col global dofs of the elem
                    ,elemgdof_col_end = elemgdofs->end() 				// same
                                ;elemgdof_col < elemgdof_col_end 		// loop over row elements from diagonal to the end
                                            ;++matloc_col
                                            ,++elemgdof_col)
                {
                    // if dof is constrained - go to next component
                    if (*elemgdof_col == Model::DOF_IS_CONSTRAINED) continue;  	// assamble only matrix columns, which correspond not constrained global dofs
                    
                    // update matrix Slices
                    elem_smatrix.update_from(matloc_col);
                    // assemble global matrix for each harmonic
                    ModelTraits::place_element_into_matrix(matrix,elem_smatrix,
                            matrix_row, *elemgdof_col-1, assemble.ndofs);
                    
                } // loop over column of the current matrix row

            } // Loop over harmonics rows

			/* store global load vector */
            // update load Slices
            sload.update_from(load.begin() + *elemgdof_row-1); // *elemgdof_row-1 == current dof number
            elem_sload.update_from(load_dof);
            // assemble global load for each harmonic
            sload += elem_sload;
        
        } // loop over inital matrix rows
		
		
		++elem; ++elemgdofs;

	} // Loop over all elements
	
}



/* Assemble extendend Jacobi as Eigen::SparseMatrix */
void ModelTraits::assemble(const Model& model, const Assemble& assemble
                            ,       Eigen::SparseMatrix<double>&   matrix      /* global system Jacobi matrix */
                            ,       math::vector_t<double,1>&   load        /* global internal system load vector */
                            , const math::vector_t<double,1>&   u           /* time domain displacement */
                            , const math::vector_t<double,1>&   dudt        /* time domain velocity */
                            , const math::vector_t<double,1>&   d2udt2      /* time domain acceleration */
                            , const math::vector_const_slice<double>& q           /* frequency domain displacement */
                            , const ::npath::DFT&               dft         /* DFT transformer */
                            , double                            freq        /* current frequency */
                            ,       math::vector_t<double,3>&   buffer_dft_matrix
                            ,       math::vector_t<double,1>&   buffer_dft_vector
                            , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                    math::vector_t<double,2>&   /* element local matrix */
                                            ,      math::vector_t<double,1>&    /* element local load */
                                            ,      math::vector_t<double,1>&    /* element extendent matrix column */ 
                                            ,const math::vector_t<double,1>&    /* element property */ 
                                            ,const math::vector_t<double,1>&    /* element material */
                                            ,const math::vector_t<double,1>&    /* element u */
                                            ,const math::vector_t<double,1>&    /* element dudt */
                                            ,const math::vector_t<double,1>&    /* element d2udt2 */
                                            ,const math::vector_const_slice<double>& /* element q */  
                                            ,const ::npath::DFT&                /* dft */ 
                                            ,double                             /* freq */ 
                                            ,      math::vector_t<double,3>&    /* buffer_dft_matrix */
                                            ,      math::vector_t<double,1>&    /* buffer_dft_vector */) const
                            )
{
	/*Iterators*/	
	decltype(assemble.elem_load.begin())  		load_dof;								/* Element load vector */
    decltype(assemble.elem_load2.begin())  		load2_dof,load2_hdof;					/* Element load2 vector current dof and corresponding harmonic */
    
 
	decltype(assemble.elem_matrix.begin()) 		matloc_row, matloc_hrow;				/* Element matrix row values: inital (zeros harmonic) and current harmonic */
	decltype(assemble.elem_matrix[0].begin()) 	matloc_col;				                /* Element matrix col values */
    
	auto 										elem = model.elements.begin()			/* Element*/
											  , elem_end = model.elements.end();	
	
	auto 										elemgdofs = assemble.elems_dofs.begin();/* Element dofs */
	decltype(assemble.elems_dofs[0].begin()) 	elemgdof_row, elemgdof_row_end 			/* Element dof for loop over element matrix rows */
											  , elemgdof_col, elemgdof_col_end;  		/*             for loop over element matrix cols */
	
    /*Slices. In time domain assemble needs loop over dofs,
    in frequency domain - loop over harmonics for each dofs is needed also.
    `s` - means Slice 
    This slices containts values in frequency domain of certain dof and all harmonics*/
    // const slices step, equal to global ndofs
    size_t hndofs = dft.frequency_size();
    math::Slice sload(  load.begin(), load.begin()+hndofs, assemble.ndofs);
    
    // step equal to element ndofs
    math::vector_slice<double> elem_sload;
    math::vector_slice<double> elem_smatrix;
    /* Other variables */
    size_t elem_ndofs, elem_hndofs, matrix_row;

	/* Loop over all elements */
	while (elem != elem_end) {
        elem_ndofs = (*elem)->ndofs();
        elem_hndofs = elem_ndofs * dft.frequency_basic_size();
		/* store element state vectors */
        ModelTraits::store_element_state_vectors(assemble.elem_displ.begin()	// where store to
												,assemble.elem_vel.begin()
												,assemble.elem_accel.begin()
                                                ,assemble.elem_displ2.begin()
												,elemgdofs->begin()				// with dofs id
												,elemgdofs->end()
												,u								// store from there
												,dudt
												,d2udt2
                                                ,q
                                                ,dft
                                                ,elem_ndofs);
        math::Slice elem_hdispl(assemble.elem_displ2.begin(),
                                assemble.elem_displ2.begin()+elem_hndofs);
        /* calculate element matrix and load vector in frequency domain */
		((*elem)->*element_matrix_load)( assemble.elem_matrix 					// where to store matrix
										,assemble.elem_load						// where to store load vector
                                        ,assemble.elem_load2
										,model.properties[(*elem)->propID]		// element property
										,model.materials[ (*elem)->matlID]		// element material
										,assemble.elem_displ					// element displacement vector
										,assemble.elem_vel
										,assemble.elem_accel
                                        ,elem_hdispl
                                        ,dft
                                        ,freq
                                        ,buffer_dft_matrix
                                        ,buffer_dft_vector);
#if 0
        if ((*elem)->nodes_dofs == 0) {
            // std::cout << "1st elem: J = \n" << assemble.elem_matrix << '\n';
            // << "\ndr/dw = " << assemble.elem_load2 << '\n';

            /* nnumerical tangent matrix calculation */
            auto num_matrix = math::zeros<double>(assemble.elem_matrix);
            auto num_load = math::zeros<double>(assemble.elem_load);
            auto num_load2 = math::zeros<double>(assemble.elem_load);

            dft.time_domain(assemble.elem_displ2,freq,assemble.elem_displ,assemble.elem_vel,assemble.elem_accel,elem_ndofs);
            (*elem)->frequency_Load( num_load						// where to store load vector
                                    ,model.properties[(*elem)->propID]		// element property
                                    ,model.materials[ (*elem)->matlID]		// element material	
                                    ,assemble.elem_displ					// element displacement vector
                                    ,assemble.elem_vel
                                    ,assemble.elem_accel
                                    ,dft
                                    ,freq
                                    ,buffer_dft_vector);
            double dx = 1e-5;
            size_t dof;
            for (dof = 0; dof < assemble.elem_displ2.size(); ++dof) {
                assemble.elem_displ2[dof] += dx;
                dft.time_domain(assemble.elem_displ2,freq,assemble.elem_displ,assemble.elem_vel,assemble.elem_accel,elem_ndofs);
                math::fill(num_load2,0.0);
                (*elem)->frequency_Load( num_load2						// where to store load vector
                                    ,model.properties[(*elem)->propID]		// element property
                                    ,model.materials[ (*elem)->matlID]		// element material	
                                    ,assemble.elem_displ					// element displacement vector
                                    ,assemble.elem_vel
                                    ,assemble.elem_accel
                                    ,dft
                                    ,freq
                                    ,buffer_dft_vector);
                for (size_t i = 0; i < assemble.elem_displ2.size(); ++i) {
                    num_matrix[i][dof] = (num_load2[i]-num_load[i])/dx;
                }
                assemble.elem_displ2[dof] -= dx;
            }
            // dr/dw
            dft.time_domain(assemble.elem_displ2,freq+dx,assemble.elem_displ,assemble.elem_vel,assemble.elem_accel,elem_ndofs);
            math::fill(num_load2,0.0);
            (*elem)->frequency_Load( num_load2						// where to store load vector
                                    ,model.properties[(*elem)->propID]		// element property
                                    ,model.materials[ (*elem)->matlID]		// element material	
                                    ,assemble.elem_displ					// element displacement vector
                                    ,assemble.elem_vel
                                    ,assemble.elem_accel
                                    ,dft
                                    ,freq+dx
                                    ,buffer_dft_vector);
            auto num_drdw = (num_load2-num_load)/dx;

            // std::cout << "num J = \n" << num_matrix << '\n'
            
            std::cout << "# |dr/dw - num dr/dw| = " << math::norm(assemble.elem_load2 - num_drdw)/math::norm(num_drdw)
            << ", |dr/dw| = " << math::norm(assemble.elem_load2)
            << ", |num dr/dw| = " << math::norm(num_drdw) << std::endl;
            // << "\nnum dr/dw = " << num_drdw << std::endl;
            assemble.elem_load2 = num_drdw;
            assemble.elem_matrix = num_matrix;
        }
#endif
		/* store element matrix and load vector to global matrix and load vector in frequency domain 
         like it in time domain */
        // update element matrix and load slices
        // - because in inner loops this slices remain const step,
        // equal to element ndofs
        elem_sload.new_slice(assemble.elem_load.begin()
                            ,assemble.elem_load.begin()
                                 +elem_ndofs*dft.frequency_basic_size()
                            ,elem_ndofs);
        elem_smatrix.new_slice(assemble.elem_matrix[0].begin(),
                               assemble.elem_matrix[0].begin()
                                   +elem_ndofs*dft.frequency_basic_size(),
                               elem_ndofs);
		/* loop over inital matrix rows (like in time domain) */
		for (matloc_row         = assemble.elem_matrix.begin()    	// iterator on elem matrix row
		    ,elemgdof_row       = elemgdofs->begin()   				// iterator on global dofs, corresponding row global dofs of the elem
			,elemgdof_row_end   = elemgdofs->end() 					// same
			,load_dof           = assemble.elem_load.begin()        // iterator on elem load vector
            ,load2_dof          = assemble.elem_load2.begin()       // iterator on elem load2 vector (here load2 is extendent matrix column)
		   				;elemgdof_row < elemgdof_row_end 			// Loop over all matrix rows
									;++matloc_row
									,++elemgdof_row
									,++load_dof
                                    ,++load2_dof)
		{

			// if dof is constrained - go to next matrix row and vector component
			if (*elemgdof_row == Model::DOF_IS_CONSTRAINED) continue;     					// assamble only matrix rows, which correspond not constrained global dofs
            
            /* Loop over harmonics rows */
            for (size_t h = 0; h < dft.frequency_basic_size(); ++h) {
                matloc_hrow = matloc_row + h*elem_ndofs;
                matrix_row = *elemgdof_row-1 + h*assemble.ndofs;

                load2_hdof = load2_dof + h*elem_ndofs;

                // update matrix Slices end bounds
                elem_smatrix.update_to(matloc_hrow->end());
                // loop over column of the current matrix row
                // loop over all element - assume that matrix is not symmetric
                for (matloc_col = matloc_hrow->begin()  					// iterator on column element
                    ,elemgdof_col = elemgdofs->begin()					// iterator on global dofs, corresponding col global dofs of the elem
                    ,elemgdof_col_end = elemgdofs->end() 				// same
                                ;elemgdof_col < elemgdof_col_end 		// loop over row elements from diagonal to the end
                                            ;++matloc_col
                                            ,++elemgdof_col)
                {
                    // if dof is constrained - go to next component
                    if (*elemgdof_col == Model::DOF_IS_CONSTRAINED) continue;  	// assamble only matrix columns, which correspond not constrained global dofs
                    
                    // update matrix Slices
                    elem_smatrix.update_from(matloc_col);
                    // assemble global matrix for each harmonic
                    ModelTraits::place_element_into_matrix(matrix,elem_smatrix,
                            matrix_row, *elemgdof_col-1, assemble.ndofs);

                
                } // loop over column of the current matrix row
                
                /* store extendent column */
                matrix.coeffRef(matrix_row, hndofs) += *load2_hdof;

            } // Loop over harmonics rows

			/* store global load vector */
            // update load Slices
            sload.update_from(load.begin() + *elemgdof_row-1); // *elemgdof_row-1 == current dof number
            elem_sload.update_from(load_dof);
            // assemble global load for each harmonic
            sload += elem_sload;
        
        } // loop over inital matrix rows
		// if ((*elem)->nodes_dofs == 0) {
        //     std::cout << "\n...1st elem: Jglobal = \n" << matrix << '\n';
        // }
		
		++elem; ++elemgdofs;

	} // Loop over all elements
	
} // Assemble extendend Jacobi as Eigen::SparseMatrix


size_t Assemble::get_max_elem_dofs() {
    return max_elem_dofs;
}


void Assemble::elem_matrix_size(size_t sz) {
    elem_matrix = math::zeros<double>(sz,sz);
    elem_matrix2 = math::zeros<double>(sz,sz);
}

void Assemble::elem_load_size(size_t sz) {
    elem_load = math::zeros<double>(sz);
}

void Assemble::elem_state_size(size_t sz) {
    elem_displ  = math::zeros<double>(sz);
    elem_vel    = math::zeros<double>(sz);
    elem_accel  = math::zeros<double>(sz);
}

void Assemble::frequency_analyses_elem_state_size(size_t freq_size,
                                                  size_t time_size) {
    size_t elem_sz = get_max_elem_dofs();
    
    elem_matrix = math::zeros<double>(freq_size*elem_sz,freq_size*elem_sz);
    elem_load = math::zeros<double>(freq_size*elem_sz);

    elem_load2 = math::zeros<double>(freq_size*elem_sz);

    elem_displ  = math::zeros<double>(time_size*elem_sz);
    elem_vel    = math::zeros<double>(time_size*elem_sz);
    elem_accel  = math::zeros<double>(time_size*elem_sz);

    elem_displ2 = math::zeros<double>(freq_size*elem_sz);
}


void ModelTraits::place_element_into_matrix(Eigen::SparseMatrix<double>& matrix
                                , const math::vector_slice<double>& local_row
                                , size_t row, size_t col, size_t step) {
    auto local = local_row.begin()
        ,local_end = local_row.end();
    while (local < local_end) {
        matrix.coeffRef(row,col) += *local;
        ++local;
        col += step;
    }
}

} // namespace fem