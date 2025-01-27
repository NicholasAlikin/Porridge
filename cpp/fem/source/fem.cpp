#include "fem.hpp"

namespace fem {

// size_t Point::EmptyID = 0;
// size_t BaseElement::EmptyID = 0;








// /* Count all releases in nodes */
// math::vector<NodeAllReleases> releases_precomputing(math::vector<ElemReleases>& elems_releases
// 											, const math::vector<BaseElement*>& elements) {
	
	
// 	std::unordered_map<size_t,NodeAllReleases> node_releases;
// 	auto elem_release = elems_releases.begin();
// 	auto elem_release_end = elems_releases.end();
// 	decltype(elems_releases[0].nodes.begin()) node_release, node_release_end;
// 	decltype(NodeAllReleases::dofs)* gnode_dofs;
// 	decltype(gnode_dofs->begin()) gnode_dof;
// 	decltype(elems_releases[0].nodes[0].dofs.begin()) node_dof, node_dof_end;

// 	size_t node_global_id;
// 	// loop over all elements with released nodes
// 	while (elem_release != elem_release_end) {
// 		node_release = elem_release->nodes.begin();
// 		node_release_end = elem_release->nodes.end();
		
		
// 		// loop over all element released nodes
// 		while (node_release != node_release_end) {
// 			node_global_id = elements[elem_release->elem_id]->nodes[node_release->node_localid];
// 			gnode_dofs = &node_releases[node_global_id].dofs;

// 			node_release->shift = node_releases[node_global_id].dofs.size()/BaseNode::DOFS;

// 			gnode_dofs->resize((node_release->shift+1)*BaseNode::DOFS);

// 			gnode_dof = gnode_dofs->begin() + node_release->shift*BaseNode::DOFS;
// 			node_dof = node_release->dofs.begin();
// 			node_dof_end = node_release->dofs.end();

// 			while (node_dof != node_dof_end) {
// 				*gnode_dof = *node_dof;
// 				++node_dof; ++gnode_dof;
// 			}
// 			++node_release;
// 		}
// 		++elem_release;
// 	}

// 	math::vector<NodeAllReleases> node_releases_vec(node_releases.size());
	
// 	auto node = node_releases.begin(), node_end = node_releases.end();
// 	auto node_vec = node_releases_vec.begin();
// 	while (node != node_end) {
// 		node->second.node_globalid = node->first;
// 		*node_vec = std::move(node->second);
// 		++node; ++node_vec;
// 	}

// 	std::sort(node_releases_vec.begin(),node_releases_vec.end());
// 	std::sort(elems_releases.begin(), elems_releases.end());
	
// 	return node_releases_vec;
// }

// void Model::parse_nodal_data(const math::vector<Node>& nodes_info
//                           , const math::vector<NodeConstraint>& constraints_info
// 						  , math::vector<BaseElement*>& elements
// 						  , const math::vector<NodeAllReleases>& releases
// 						  , const math::vector<ElemReleases>& elems_releases)
// {
//     // matrix with nodes global dofs
//     math::vector_t<size_t,2> nodes_dofs = math::zeros<size_t>(nodes_info.size(),BaseNode::DOFS);
    
//     // iterators
//     auto gnode = nodes_dofs.begin(), node_end = nodes_dofs.end();
// 	size_t gnodeID = 0;
//     typename math::vector<size_t>::iterator gdof, gdof_end;

//     auto node_info = nodes_info.begin();
//     auto constraint_node = constraints_info.begin();
// 	typename decltype(NodeConstraint::dofs)::const_iterator constraint_dof;
    
// 	auto release = releases.begin()
// 		,release_end = releases.end();
// 	decltype(releases[0].dofs.begin()) released_dof,released_dof_end;
// 	// set nodes_dofs matrix by global dofs if current dof is not constrained
// 	// constrained dofs value = 0
//     size_t global_id = 1;
//     while (gnode != node_end) {
//         gdof = gnode->begin();
// 		gdof_end = gnode->end();
        
// 		// if dof is constained, gdof=0
// 		if (node_info->id == constraint_node->id) {
// 			constraint_dof = constraint_node->dofs.begin();
// 			while (gdof != gdof_end) {
// 				if (*constraint_dof == 0) {
// 					*gdof = global_id;
// 					++global_id;
// 				}
// 				++gdof; ++constraint_dof;
// 			}
// 			++constraint_node;
//         } else {
// 			while (gdof != gdof_end) {
// 				*gdof = global_id;
// 				++global_id; ++gdof;
// 			}
//         }

// 		// Released dof
// 		if ((release < release_end) && (gnodeID == release->node_globalid)) {
// 			// resize vector, iterator gnode not invalidate
// 			gnode->resize(BaseNode::DOFS + release->dofs.size());
			
// 			gdof = gnode->begin()+BaseNode::DOFS; // skip already desined dofs
// 			released_dof = release->dofs.begin();
// 			released_dof_end = release->dofs.end();
// 			// loop over released dof
// 			while (released_dof != released_dof_end) {
				
// 				if (*released_dof == NodeAllReleases::IsNotReleased) {
// 					*gdof = 0;
// 				} else {
// 					*gdof = global_id;
// 					++global_id;
// 				}
// 				++released_dof; ++gdof;
// 			}
			
// 			++release;
// 		}


// 		++gnode; ++node_info; ++gnodeID;
//     }
	
// 	// Calculate elements properties
// 	auto elem = elements.begin(), elem_end = elements.end();
// #ifdef BEAM_RELEASES
// 	auto elem_release = elems_releases.begin()
// 	    ,elem_release_end = elems_releases.end();
// #endif

// 	while (elem != elem_end) {
		
// 		#ifdef BEAM_RELEASES
// 		if ((elem_release < elem_release_end) && ((*elem)->nodes_dofs == elem_release->elem_id)) {
// 			(*elem)->calc_parameters(nodes_info,*elem_release);
// 			++elem_release;
			
// 			++elem; // don`t forget to increment element iterator
// 			continue;
// 		}
// 		#endif
		
// 		(*elem)->calc_parameters(nodes_info);
		
// 		++elem;
// 	}

// 	// Added dofs according to beam releases dofs
// 	// new dofs are added to back of the node already defined dofs
// 	// node_dofs = [g1,g2,g3,g4,g5,g6, r1,r2,r3,r4,r5,r6]
// 	// where gi - defined dofs
// 	// ri = 0 if dof is not released
// 	// ri = gdofID if dof is released
// 	// so element reliased dofs = gi + BaseNode::DOFS = gi + 6
	


// 	Assemble GDofs(global_id-1,std::move(nodes_dofs)); // return value
// 	return GDofs;								 // optimization
// }






// math::vector_t<double,2> assemble(Assemble& GDofs, const math::vector<BaseElement*>& elements
// 				, const math::vector_t<double,2>& properties
// 				, const math::vector_t<double,2>& materials
// 				, math::vector_t<double>& band
// 				, math::vector_t<double,2> (BaseElement::* element_matrix)(const math::vector<double>&,const math::vector<double>&) const)
// {
// 	// 2. fill global matrix
// 	// TODO
// 	math::vector_t<double,2> matrix = math::zeros<double>(GDofs.ndofs,GDofs.ndofs);
// 	// band = math::zeros<double>(math::sum(GDofs.colhs));
	
// 	math::vector_t<double,2> matrix_local; // element local matrix
// 	decltype(matrix_local.begin()) matloc_row,matloc_row_end;
// 	decltype(matrix_local[0].begin()) matloc_col,matloc_col_end;

// 	auto elem = elements.begin(), elem_end = elements.end();
// 	auto elemgdofs = GDofs.elems_dofs.begin(); // iterator on current element dofs
// 	decltype(GDofs.elems_dofs[0].begin()) elemgdof_row, elemgdof_row_end
// 									, elemgdof_col, elemgdof_col_end; // iterator on current element dof
	

// 	decltype(elements[0]->nodes.begin()) elemnode, elemnode_end; // iterator on element node
// 	decltype(GDofs.nodes_dofs[0].begin()) nodegdof, nodegdof_end; // iterator on element global dof
// 	size_t matloc_row_num;
// 	// loop over elements
// 	while (elem != elem_end) {
// 		// calculate element matrix (stiffness, mass)
// 		matrix_local = ((*elem)->*element_matrix)(properties[(*elem)->propID]
// 										         , materials[(*elem)->matlID]);
		
// 		// loop over element matrix
// 		// loop over matrix rows
// 		for (matloc_row = matrix_local.begin()    // iterator on elem matrix row
// 		    ,matloc_row_end = matrix_local.end()  // same
// 		    ,elemgdof_row = elemgdofs->begin()    // iterator on global dofs, corresponding row global dofs of the elem
// 			,matloc_row_num = 0					  // elem matrix count - to go throw the upper triangular part of the matrix only
// 					;matloc_row != matloc_row_end // loop over all matrix rows
// 									;++matloc_row
// 									,++elemgdof_row
// 									,++matloc_row_num) {
// 			if (*elemgdof_row == 0) continue;     // assamble only matrix rows, which correspond not constrained global dofs
// 			// loop over columnt of current row
// 			for (matloc_col = matloc_row->begin() + matloc_row_num  // iterator on column element - only upper triangular part
// 			    ,matloc_col_end = matloc_row->end()                 // same
// 				,elemgdof_col = elemgdofs->begin() + matloc_row_num // iterator on global dofs, corresponding col global dofs of the elem
// 						;matloc_col != matloc_col_end // loop over row elements from diagonal to the end
// 								;++matloc_col
// 								,++elemgdof_col) {
// 				if (*elemgdof_col == 0) continue;  // assamble only matrix columns, which correspond not constrained global dofs
				
// 				// global matrix like 2d array - stored only upper triangular part, lower triangular is symmetric
// 				matrix[*elemgdof_row-1][*elemgdof_col-1] += *matloc_col;
				
// 				// global matrix like 1d array - stored only upper triangular part
// 				// GDofs.diags - addresses of the matrix diagonal elements
// 				// &Kij = &diag + col_num - row_num
// 				band[(GDofs.diags[*elemgdof_col-1]
// 				     + *elemgdof_col) - *elemgdof_row] += *matloc_col;
				
// 			}
// 		}
		
// 		++elem; ++elemgdofs;
// 	}

// 	return matrix;
// }


// void assembleNL(const Assemble& GDofs, const math::vector<BaseElement*>& elements
// 				, const math::vector_t<double,2>& properties
// 				, const math::vector_t<double,2>& materials
// 				, math::vector<double>& band
// 				, math::vector<double>& load
// 				, void (BaseElement::* element_matrix_load)(math::vector_t<double,2>& /* K */
// 													 ,      math::vector_t<double,1>& /* internal_load */ 
// 													 ,const math::vector_t<double,1>& /* property */ 
// 													 ,const math::vector_t<double,1>& /* material */ 
// 													 ,const math::vector_t<double,3>& /* Rsum */ 
// 													 ,const math::vector_t<double,1>& /* q */  ) const
// 				, const math::vector<double>& q
// 				, const math::vector_t<double,3>& Rsum
// 				, math::vector_t<double,2>& matrix_local
// 				, math::vector<double>& load_local
// 				, math::vector<double>& q_elem)
// {
// 	// 2. fill global matrix
// 	/*
	



// 	* q - 
// 		Global displacement vector.
// 	* Rsum -
// 		Array of total rotation tensors of nodes, stored by nodes global IDs.

// 	*/
// 	for (double& x: band)
// 		x = 0.0;
// 	for (double& x: load)
// 		x = 0.0;
		
// 	// TODO
// 	// math::vector_t<double,2> matrix = math::zeros<double>(GDofs.ndofs,GDofs.ndofs);
// 	// band = math::zeros<double>(math::sum(GDofs.colhs));
	
// 	// math::vector_t<double,2> matrix_local; // element local matrix
// 	// math::vector<double> load_local, // element local load: internal
// 	// 						q_elem;  // element displacements
// 	decltype(q_elem.begin()) q_elem_dof,q_elem_dof_end;
// 	decltype(load_local.begin()) load_dof;

// 	decltype(matrix_local.begin()) matloc_row,matloc_row_end;
// 	decltype(matrix_local[0].begin()) matloc_col,matloc_col_end;

// 	auto elem = elements.begin(), elem_end = elements.end();
// 	auto elemgdofs = GDofs.elems_dofs.begin(); // iterator on current element dofs
// 	decltype(GDofs.elems_dofs[0].begin()) elemgdof_row, elemgdof_row_end
// 									, elemgdof_col, elemgdof_col_end; // iterator on current element dof
	

// 	decltype(elements[0]->nodes.begin()) elemnode, elemnode_end; // iterator on element node
// 	decltype(GDofs.nodes_dofs[0].begin()) nodegdof, nodegdof_end; // iterator on element global dof
// 	size_t matloc_row_num;
// 	// loop over elements
// 	while (elem != elem_end) {

// 		// store element displacements
// 		// loop over element dofs
// 		for (q_elem_dof = q_elem.begin()     		// iterator on elem displacements
// 			,elemgdof_row = elemgdofs->begin()		// iterator in global dofs nodes_dofs, corresponding elem global dofs
// 			,elemgdof_row_end = elemgdofs->end()
// 					;elemgdof_row < elemgdof_row_end   //
// 							;++q_elem_dof
// 							,++elemgdof_row) {
// 			if (*elemgdof_row == 0) {
// 				*q_elem_dof = 0.0;
// 				continue;       // if dof is constained its displacement = 0
// 			}
// 			*q_elem_dof = q[*elemgdof_row-1];		// store element displacement value of *elemgdof_row dof nodes_dofs
// 		}

// 		// calculate element matrix (stiffness, mass) and load (internal, inertion)
// 		((*elem)->*element_matrix_load)( matrix_local 
// 										,load_local
// 										,properties[(*elem)->propID]
// 										,materials[(*elem)->matlID]
// 										,Rsum
// 										,q_elem);
// 		// std::cout << "matrix_local = \n" << matrix_local << std::endl;
// 		// loop over element matrix
// 		// loop over matrix rows
// 		for (matloc_row = matrix_local.begin()    // iterator on elem matrix row
// 		    ,elemgdof_row = elemgdofs->begin()    // iterator on global dofs, corresponding row global dofs of the elem
// 			,elemgdof_row_end = elemgdofs->end()    // same
// 			,matloc_row_num = 0					  // elem matrix count - to go throw the upper triangular part of the matrix only
// 			,load_dof = load_local.begin()			  // iterator on elem load vector
// 		   				;elemgdof_row < elemgdof_row_end // loop over all matrix rows
// 									;++matloc_row
// 									,++elemgdof_row
// 									,++matloc_row_num
// 									,++load_dof) {
// 			if (*elemgdof_row == 0) continue;     // assamble only matrix rows, which correspond not constrained global dofs

// 			// loop over columnt of current row
// 			for (matloc_col = matloc_row->begin() + matloc_row_num  // iterator on column element - only upper triangular part
// 			    ,elemgdof_col = elemgdofs->begin() + matloc_row_num // iterator on global dofs, corresponding col global dofs of the elem
// 				,elemgdof_col_end = elemgdofs->end() // same
// 							;elemgdof_col < elemgdof_col_end // loop over row elements from diagonal to the end
// 										;++matloc_col
// 										,++elemgdof_col) {
// 				if (*elemgdof_col == 0) continue;  // assamble only matrix columns, which correspond not constrained global dofs


// 				// global matrix like 2d array - stored only upper triangular part, lower triangular is symmetric
// 				// matrix[*elemgdof_row-1][*elemgdof_col-1] += *matloc_col;
// 				// matrix[*elemgdof_col-1][*elemgdof_row-1] += *matloc_col;
				
// 				// global matrix like 1d array - stored only upper triangular part
// 				// GDofs.diags - addresses of the matrix diagonal elements
// 				// &Kij = &diag + col_num - row_num
// 				band[(GDofs.diags[*elemgdof_col-1]
// 				     + *elemgdof_col) - *elemgdof_row] += *matloc_col;
				
// 			}

// 			// store global load vector
// 			load[*elemgdof_row-1] += *load_dof;
// 		}
		
		
// 		++elem; ++elemgdofs;
// 	}
	
// 	// return matrix;
// }

// void assembleNL(const Assemble& GDofs, const math::vector<BaseElement*>& elements
// 				, const math::vector_t<double,2>& properties
// 				, const math::vector_t<double,2>& materials
// 				, math::vector_t<double,2>& matrix
// 				, math::vector<double>& load
// 				, void (BaseElement::* element_matrix_load)(math::vector_t<double,2>& /* K */
// 													 ,      math::vector_t<double,1>& /* internal_load */ 
// 													 ,const math::vector_t<double,1>& /* property */ 
// 													 ,const math::vector_t<double,1>& /* material */ 
// 													 ,const math::vector_t<double,3>& /* Rsum */ 
// 													 ,const math::vector_t<double,1>& /* q */  ) const
// 				, const math::vector<double>& q
// 				, const math::vector_t<double,3>& Rsum)
// {
// 	// 2. fill global matrix
// 	/*
	



// 	* q - 
// 		Global displacement vector.
// 	* Rsum -
// 		Array of total rotation tensors of nodes, stored by nodes global IDs.

// 	*/

	

// 	for (auto& xr: matrix) {
// 		for (double& xc: xr) {
// 			xc = 0.0;
// 		}
// 	}
// 	// matrix = math::zeros<double>(GDofs.ndofs, GDofs.ndofs);
// 	for (double& x: load)
// 		x = 0.0;
		
	
// 	// TODO
// 	// math::vector_t<double,2> matrix = math::zeros<double>(GDofs.ndofs,GDofs.ndofs);
// 	// band = math::zeros<double>(math::sum(GDofs.colhs));
	
// 	math::vector_t<double,2> matrix_local; // element local matrix
// 	math::vector<double> load_local, // element local load: internal
// 							q_elem;  // element displacements
// 	decltype(q_elem.begin()) q_elem_dof,q_elem_dof_end;
// 	decltype(load_local.begin()) load_dof;

// 	decltype(matrix_local.begin()) matloc_row,matloc_row_end;
// 	decltype(matrix_local[0].begin()) matloc_col,matloc_col_end;

// 	auto elem = elements.begin(), elem_end = elements.end();
// 	auto elemgdofs = GDofs.elems_dofs.begin(); // iterator on current element dofs
// 	////
	
// 	////
// 	decltype(GDofs.elems_dofs[0].begin()) elemgdof_row, elemgdof_row_end
// 									, elemgdof_col, elemgdof_col_end; // iterator on current element dof
	

// 	decltype(elements[0]->nodes.begin()) elemnode, elemnode_end; // iterator on element node
// 	decltype(GDofs.nodes_dofs[0].begin()) nodegdof, nodegdof_end; // iterator on element global dof
// 	size_t matloc_row_num;
// 	// loop over elements
// 	while (elem < elem_end) {

// 		// store element displacements
// 		q_elem = math::zeros<double>((*elem)->ndofs());
// 		// loop over element dofs
		
		
		
// 		for (q_elem_dof = q_elem.begin()		// iterator on elem displacements
// 			,q_elem_dof_end = q_elem.end()
// 		    ,elemgdof_row = elemgdofs->begin()		// iterator in global dofs nodes_dofs, corresponding elem global dofs
// 					;q_elem_dof != q_elem_dof_end   //
// 							;++q_elem_dof
// 							,++elemgdof_row)
// 		{
// 			if (*elemgdof_row == 0) continue;       // if dof is constained its displacement = 0
// 			*q_elem_dof = q[*elemgdof_row-1];		// store element displacement value of *elemgdof_row dof nodes_dofs
// 		}
// 		// calculate element matrix (stiffness, mass) and load (internal, inertion)
// 		((*elem)->*element_matrix_load)( matrix_local 
// 										,load_local
// 										,properties[(*elem)->propID]
// 										,materials[(*elem)->matlID]
// 										,Rsum
// 										,q_elem);
		
// 		// loop over element matrix
// 		// loop over matrix rows
// 		for (matloc_row = matrix_local.begin()    // iterator on elem matrix row
// 		    ,matloc_row_end = matrix_local.end()  // same
// 		    ,elemgdof_row = elemgdofs->begin()    // iterator on global dofs, corresponding row global dofs of the elem
// 			,matloc_row_num = 0					  // elem matrix count - to go throw the upper triangular part of the matrix only
// 			,load_dof = load_local.begin()			  // iterator on elem load vector
// 		   				;matloc_row != matloc_row_end // loop over all matrix rows
// 									;++matloc_row
// 									,++elemgdof_row
// 									,++matloc_row_num
// 									,++load_dof) {
// 			if (*elemgdof_row == 0) continue;     // assamble only matrix rows, which correspond not constrained global dofs

// 			// loop over columnt of current row
// 			for (matloc_col = matloc_row->begin() + matloc_row_num  // iterator on column element - only upper triangular part
// 			    ,matloc_col_end = matloc_row->end()                 // same
// 				,elemgdof_col = elemgdofs->begin() + matloc_row_num // iterator on global dofs, corresponding col global dofs of the elem
// 							;matloc_col != matloc_col_end // loop over row elements from diagonal to the end
// 										;++matloc_col
// 										,++elemgdof_col) {
// 				if (*elemgdof_col == 0) continue;  // assamble only matrix columns, which correspond not constrained global dofs


// 				// global matrix like 2d array - stored only upper triangular part, lower triangular is symmetric
// 				matrix[*elemgdof_row-1][*elemgdof_col-1] += *matloc_col;
// 				if (*elemgdof_col != *elemgdof_row)
// 					matrix[*elemgdof_col-1][*elemgdof_row-1] += *matloc_col;
				
// 				// global matrix like 1d array - stored only upper triangular part
// 				// GDofs.diags - addresses of the matrix diagonal elements
// 				// &Kij = &diag + col_num - row_num
// 				// band[(GDofs.diags[*elemgdof_col-1]
// 				//      + *elemgdof_col) - *elemgdof_row] += *matloc_col;
				
// 			}

// 			// store global load vector
// 			load[*elemgdof_row-1] += *load_dof;
// 		}
		
		
// 		++elem; ++elemgdofs;
// 	}
	
// 	// return matrix;
// }


} // namespace fem