#include "analyzes.hpp"

namespace fem {

void AnalysisTraits::setup_analysis(const Model& model, Assemble& assemble) {
    ModelTraits::parse_nodal_data(model,assemble);
    ModelTraits::assemble_precomputing(model,assemble);
}

math::vector_t<double,1> AnalysisTraits::staticLD(Model& model
								, Assemble& assemble
				                , size_t load_steps
								, double epsq
								, double epsload)
{
	/*Static analysis with large displacements, linear elasticity.
	
	1. Inital iteration like for linear static analysis
	2. Load step
	2.1. Newton-Raphson iterations while displacement and/or load convergence.

	*/

	
    // assemble_precomputing(GDofs,elements,elems_releases);
    AnalysisTraits::setup_analysis(model,assemble);
	// initialazing stiffness matrix, load vector and displacement vector
	math::vector_t<double> stif(math::sum(assemble.colhs))
                        ,  loadExt(assemble.ndofs)
                        ,  q(assemble.ndofs)
                        ,  loadInt(assemble.ndofs)
                        ,  dq(assemble.ndofs)
						,  dload;
    /*Temp variables*/
	math::vector<double> temp_theta(3);
	math::vector_t<double,2> temp_rotTensor = math::zeros<double>(3,3);
	math::vector_t<double,2> temp_Rsumi = math::zeros<double>(3,3);
	
    // Calculate exact load vector
    // store_load_vector(loadExt,GDofs,loads_info);
	ModelTraits::store_load_vector(model,assemble,loadExt);
    double loadExt_abs = norm(loadExt);
	// std::cout << "GDofs.nodes_dofs = \n" << GDofs.nodes_dofs << std::endl;
	// std::cout << "GDofs.ndofs = " << GDofs.ndofs << std::endl;
	// std::cout << "loadExt = \n" << loadExt << std::endl;
    
	
    // Inital value of q equal to 0, so zeros iteration like static linear analysis.
	double cur_loadExt_abs = loadExt_abs/load_steps;
    loadExt *= cur_loadExt_abs/loadExt_abs;
    
	math::vector_t<double,3> Rsum = AnalysisTraits::setup_Rsum(model);

	// std::cout << "Rsum = \n" << Rsum << std::endl;
	// throw 1;
    
	// assemble linear like system
	// assembleNL(GDofs,elements,properties,materials,stif,loadInt,&BaseElement::tangentStiffness_innerLoad,q,Rsum
	// 				,GDofs.elem_matrix,GDofs.elem_load, GDofs.elem_displ);
	ModelTraits::do_assemble_nonlinear(model,assemble,stif,loadInt,q,Rsum,
										&BaseElement::tangentStiffness_innerLoad);
	// std::cout << "q = \n" << q << std::endl;
	// std::cout << "loadInt = \n" << loadExt << std::endl;
	// std::cout << "Ktmp = \n" << Ktmp << std::endl;
	// std::cout << "stif = \n" << stif << std::endl;

    std::cout << "#norm(loadExt) = " << math::norm(loadExt) << std::endl;
	std::cout << "#norm(loadInt) = " << math::norm(loadInt) << std::endl;
	size_t subiter = 1;
	dload = loadExt;
	auto LT = math::zeros<double>(stif);
	auto D = math::zeros<double>(assemble.ndofs);
	math::vector_t<double,2> rotTensor = math::zeros<double>(3,3);
	math::vector_t<double,2> Rtemp(rotTensor);

    math::solve_ldlt(stif,assemble.diags,loadExt,dq,LT,D,assemble.ndofs);
	q += dq;
	while ((math::norm(dq) > epsq) || (math::norm(dload) > epsload)) {
		AnalysisTraits::update_Rsum(model,assemble,Rsum,q,temp_theta,temp_rotTensor,temp_Rsumi);
		math::fill(stif.begin(),stif.end(), 0.0);
    	math::fill(loadInt.begin(),loadInt.end(), 0.0);
		ModelTraits::do_assemble_nonlinear(model,assemble,stif,loadInt,q,Rsum,&BaseElement::tangentStiffness_innerLoad);

		dload = loadExt - loadInt;
		math::solve_ldlt(stif,assemble.diags,dload,dq,LT,D,assemble.ndofs);
		q += dq;
		std::cout << "\t# subiter "<< subiter <<", |dq| = " << math::norm(dq) << ", |dload| = " << math::norm(dload) << std::endl;
		++subiter;
	}
	std::cout << "#\titer "<< 0 <<", subiters " << subiter << std::endl;

    for (size_t iter = 0; iter < load_steps-1; ++iter) {
		loadExt *= (1.0 + 1.0/(iter+1));
		std::cout << "#norm(loadExt) = " << math::norm(loadExt) << std::endl;
		
		AnalysisTraits::update_Rsum(model,assemble,Rsum,q,temp_theta,temp_rotTensor,temp_Rsumi);
		math::fill(stif.begin(),stif.end(), 0.0);
    	math::fill(loadInt.begin(),loadInt.end(), 0.0);
		ModelTraits::do_assemble_nonlinear(model,assemble,stif,loadInt,q,Rsum,&BaseElement::tangentStiffness_innerLoad);

		dload = loadExt - loadInt;
		math::solve_ldlt(stif,assemble.diags,dload,dq,LT,D,assemble.ndofs);
		q += dq;
		subiter = 1;
		std::cout << "\t# subiter "<< subiter <<", |dq| = " << math::norm(dq) << ", |dload| = " << math::norm(dload) << std::endl;
		while ((math::norm(dq) > epsq) || (math::norm(dload) > epsload)) {
	
			AnalysisTraits::update_Rsum(model,assemble,Rsum,q,temp_theta,temp_rotTensor,temp_Rsumi);
			math::fill(stif.begin(),stif.end(), 0.0);
    		math::fill(loadInt.begin(),loadInt.end(), 0.0);
			ModelTraits::do_assemble_nonlinear(model,assemble,stif,loadInt,q,Rsum,&BaseElement::tangentStiffness_innerLoad);

            dload = loadExt - loadInt;
			math::solve_ldlt(stif,assemble.diags,dload,dq,LT,D,assemble.ndofs);
            q += dq;
			++subiter;
			std::cout << "\t# subiter "<< subiter <<", |dq| = " << math::norm(dq) << ", |dload| = " << math::norm(dload) << std::endl;
        }
		std::cout << "# iter "<< iter <<", subiters " << subiter << std::endl;
    }
	std::cout << "#norm(loadExt) = " << math::norm(loadExt) << std::endl;

	return q;
}


math::vector_t<double,1> AnalysisTraits::staticLD2(
											Model& model
											,Assemble& assemble
											,size_t load_steps
                							, double epsq
											, double epsload)
{
	/*Static analysis with large displacements, linear elasticity.
	
	1. Inital iteration like for linear static analysis
	2. Load step
	2.1. Newton-Raphson iterations while displacement and/or load convergence.

	*/
	
	// precomputing
	AnalysisTraits::setup_analysis(model,assemble);
    // initialazing stiffness matrix, load vector and displacement vector
	math::vector_t<double> loadExt(assemble.ndofs)
                        ,  q(assemble.ndofs)
                        ,  loadInt(assemble.ndofs)
                        ,  dq(assemble.ndofs)
						,  dload;
    math::vector_t<double,2> stif = math::zeros<double>(assemble.ndofs,assemble.ndofs);
    // Calculate exact load vector
    ModelTraits::store_load_vector(model,assemble,loadExt);
    double loadExt_abs = norm(loadExt);
	// std::cout << "GDofs.nodes_dofs = \n" << GDofs.nodes_dofs << std::endl;
	// std::cout << "GDofs.ndofs = " << GDofs.ndofs << std::endl;
	// std::cout << "loadExt = \n" << loadExt << std::endl;
    
	
    // Inital value of q equal to 0, so zeros iteration like static linear analysis.
	double cur_loadExt_abs = loadExt_abs/load_steps;
    loadExt *= cur_loadExt_abs/loadExt_abs;
    
	math::vector_t<double,3> Rsum = AnalysisTraits::setup_Rsum(model);

	// std::cout << "Rsum = \n" << Rsum << std::endl;
	// throw 1;
    
	// assemble linear like system
	ModelTraits::do_assemble_nonlinear(model,assemble,stif,loadInt,q,Rsum,&BaseElement::tangentStiffness_innerLoad);
	// std::cout << "q = \n" << q << std::endl;
	// std::cout << "loadInt = \n" << loadExt << std::endl;
	// std::cout << "Ktmp = \n" << Ktmp << std::endl;
	// std::cout << "stif = \n" << stif << std::endl;

    std::cout << "#norm(loadExt) = " << math::norm(loadExt) << std::endl;
	std::cout << "#norm(loadInt) = " << math::norm(loadInt) << std::endl;
	size_t subiter = 1;
	dload = loadExt;
	auto L = math::zeros<double>(stif);
	auto D = math::zeros<double>(assemble.ndofs);
	math::vector_t<double,2> temp_rotTensor = math::zeros<double>(3,3);
	math::vector_t<double,2> temp_Rsumi(temp_rotTensor);
	math::vector<double> temp_theta(3);

    math::solve_ldlt(stif,loadExt,dq,L,D,assemble.ndofs);
	q += dq;
	while ((math::norm(dq) > epsq) || (math::norm(dload) > epsload)) {
		AnalysisTraits::update_Rsum(model,assemble,Rsum,q,temp_theta,temp_rotTensor,temp_Rsumi);
		math::fill(stif.begin(),stif.end(), 0.0);
    	math::fill(loadInt.begin(),loadInt.end(), 0.0);
		ModelTraits::do_assemble_nonlinear(model,assemble,stif,loadInt,q,Rsum,&BaseElement::tangentStiffness_innerLoad);
	
		dload = loadExt - loadInt;
		math::solve_ldlt(stif,dload,dq,L,D,assemble.ndofs);
		q += dq;
		std::cout << "\t# subiter "<< subiter <<", |dq| = " << math::norm(dq) << ", |dload| = " << math::norm(dload) << std::endl;
		++subiter;
	}
	std::cout << "#\titer "<< 0 <<", subiters " << subiter << std::endl;

    for (size_t iter = 0; iter < load_steps-1; ++iter) {
		loadExt *= (1.0 + 1.0/(iter+1));
		std::cout << "#norm(loadExt) = " << math::norm(loadExt) << std::endl;
		AnalysisTraits::update_Rsum(model,assemble,Rsum,q,temp_theta,temp_rotTensor,temp_Rsumi);
		math::fill(stif.begin(),stif.end(), 0.0);
    	math::fill(loadInt.begin(),loadInt.end(), 0.0);
		ModelTraits::do_assemble_nonlinear(model,assemble,stif,loadInt,q,Rsum,&BaseElement::tangentStiffness_innerLoad);
		
		dload = loadExt - loadInt;
		math::solve_ldlt(stif,dload,dq,L,D,assemble.ndofs);
		q += dq;
		subiter = 1;
		std::cout << "\t# subiter "<< subiter <<", |dq| = " << math::norm(dq) << ", |dload| = " << math::norm(dload) << std::endl;
		while ((math::norm(dq) > epsq) || (math::norm(dload) > epsload)) {
			AnalysisTraits::update_Rsum(model,assemble,Rsum,q,temp_theta,temp_rotTensor,temp_Rsumi);
			math::fill(stif.begin(),stif.end(), 0.0);
    		math::fill(loadInt.begin(),loadInt.end(), 0.0);
			ModelTraits::do_assemble_nonlinear(model,assemble,stif,loadInt,q,Rsum,&BaseElement::tangentStiffness_innerLoad);
			
			dload = loadExt - loadInt;
			math::solve_ldlt(stif,dload,dq,L,D,assemble.ndofs);
            q += dq;
			++subiter;
			std::cout << "\t# subiter "<< subiter <<", |dq| = " << math::norm(dq) << ", |dload| = " << math::norm(dload) << std::endl;
        }
		std::cout << "# iter "<< iter <<", subiters " << subiter << std::endl;
    }
	std::cout << "#norm(loadExt) = " << math::norm(loadExt) << std::endl;

	return q;
}

void update_loadExt_vector(const Assemble& GDofs
						, math::vector<fem::NodeLoad>& loads_info
						, math::vector<double>& loadExt
						, double& cur_loadExt_norm
						, double new_loadExt_norm)
{
	auto load_node     = loads_info.begin()
		,load_node_end = loads_info.end();
		
	decltype(loads_info[0].dofs.begin()) load_dof, load_dof_end;
	
	decltype(GDofs.nodes_dofs.begin()) gnode;
	decltype(GDofs.nodes_dofs[0].begin()) gdof;
	
	while (load_node != load_node_end) {
		load_dof = load_node->dofs.begin();
		load_dof_end = load_node->dofs.end();
		gdof = GDofs.nodes_dofs[load_node->id].begin();
		while (load_dof != load_dof_end) {
			if (*gdof > 0) {
                *load_dof *= new_loadExt_norm/cur_loadExt_norm;
				loadExt[*gdof-1] = *load_dof;
			}
			++load_dof; ++gdof;
		}
		++load_node;
	}
    cur_loadExt_norm = new_loadExt_norm;

}

// void update_Rsum(const Assemble& assemble
// 				,math::vector_t<double,3>& Rsum
// 				,math::vector<double>& q
// 				,math::vector<double> theta
// 				,math::vector_t<double,2>& rotTensor
// 				,math::vector_t<double,2>& Rtemp)
// {
// 	/*For each node: Rsumi = dot(L(thetai),Rsumi)
// 	  thetai_x = thetai_y = thetai_z = 0;
// 	*/
// 	// math::vector<double> theta(3);
// 	double theta_abs_2;
// 	decltype(theta.begin()) theta_dof, theta_dof_end = theta.end();

// 	auto node = assemble.nodes_dofs.begin()
// 	   , node_end = assemble.nodes_dofs.end();
// 	decltype(assemble.nodes_dofs[0].begin()) gdof, gdof_end;
// 	auto Rsum_node = Rsum.begin();

// 	// loop over all nodes
// 	while (node < node_end) {
		
// 		// num_subnodes = node->size()/BaseNode::DOFS-1;
// 		// for (subnode = 0; subnode < num_subnodes) {

// 		// }
// 		// store theta of node
// 		theta_abs_2 = 0;
// 		for (gdof = node->begin() + 3 // dofs = [ux,uy,uz,tx,ty,tz] [ux,uy,uz,tx,ty,tz] [ux,uy,uz,tx,ty,tz]
// 			,theta_dof = theta.begin()
// 					;theta_dof < theta_dof_end
// 							;++gdof
// 							,++theta_dof)
// 		{
// 			if (*gdof == 0) {
// 				*theta_dof = 0.0;
// 				continue;
// 			}
			
// 			*theta_dof = q[*gdof-1]; // 
// 			q[*gdof-1] = 0;          //
// 			theta_abs_2 += (*theta_dof)*(*theta_dof);
// 		}

// 		// Rsumi = dot(L(thetai),Rsumi)
// 		// update Rsumi only if rotation is not 0
// 		if (theta_abs_2 > 0) {
// 			// *Rsum_node = math::dot(math::rotation_tensor(theta), *Rsum_node);
// 			math::rotation_tensor(theta,rotTensor);
// 			math::make_zeros(Rtemp);
// 			math::dot(rotTensor,*Rsum_node,Rtemp);
// 			math::swap(*Rsum_node,Rtemp);
// 		}
// 		++node;	++Rsum_node;
// 	}

// }

void AnalysisTraits::update_Rsum(const Model& model
								,const Assemble& assemble
								,math::vector_t<double,3>& Rsum
								,math::vector<double>& q
								,math::vector<double>& theta
								,math::vector_t<double,2>& temp_rotTensor
								,math::vector_t<double,2>& temp_Rsumi)
{
	/*For each node: Rsumi = dot(L(thetai),Rsumi)
	  thetai_x = thetai_y = thetai_z = 0;
	*/
	// double theta_abs_2;
	// decltype(theta.begin()) theta_dof, theta_dof_end = theta.end();

	auto node = assemble.nodes_dofs.begin()
	   , node_end = assemble.nodes_dofs.end();
	size_t node_id = 0;
	decltype(assemble.nodes_dofs[0].begin()) gdof, subdof, subdof_end;
	auto Rsum_node = Rsum.begin();

	// nodes with releases
	auto node_release = model.nodes_releases.begin()
		,node_release_end = model.nodes_releases.end();
	size_t num_subnodes, subnode;
	decltype(Rsum[0].begin()) itRsumi;


	// loop over all nodes
	for (;node < node_end; 	 ++node
							,++Rsum_node
							,++node_id)
	{
		
		/* nodes with releases
		possible node gdofs = [d1,d2,d3,d4,d5,d6] [0,d7,0,0,0,d8] [0,0,d9,0,0,0]
		*/
#ifdef BEAM_RELEASES
		if ((node_release < node_release_end) && (node_id == node_release->node_globalid)) {
			// loop over only 'sub' nodes!
			num_subnodes = node->size()/BaseNode::DOFS; // '0' is 'main' node
			for (subnode = 1; subnode < num_subnodes; ++subnode) {
				/*Go thought the sub and main nodes,
				  take only not released dofs from the sub node*/
				
				AnalysisTraits::update_Rsumi_releases(q
									 ,Rsum_node->begin() + 3*subnode		// Rsumi
									 ,node->begin() + 3					// gdof
									 ,node->begin() + 3 + BaseNode::DOFS*subnode		// subdof
									 ,theta
									 ,temp_rotTensor
									 ,temp_Rsumi);

			}
			// main node - subnode = 0;
			AnalysisTraits::update_Rsumi_releases(q
								 ,Rsum_node->begin()		// Rsumi
								 ,node->begin() + 3					// gdof
								 ,theta
								 ,temp_rotTensor
								 ,temp_Rsumi);

			++node_release;
			continue;

		}
#endif

		AnalysisTraits::update_Rsumi(q
								, Rsum_node
								, node->begin() + 3
								, theta
								, temp_rotTensor
								, temp_Rsumi);
		
	}
	// std::cout << "Rsum = \n" << Rsum << std::endl;

}

void AnalysisTraits::update_Rsumi_releases(math::vector<double>& q
								, typename math::vector_t<double,2>::iterator Rsumi
								, typename math::vector<size_t>::const_iterator gdof
								, typename math::vector<size_t>::const_iterator subdof
								, 		   math::vector<double>& theta
								, math::vector_t<double,2>& temp_rotTensor
								, math::vector_t<double,2>& temp_Rsumi)
{
	auto theta_dof = theta.begin()
		,theta_dof_end = theta.end();
	double theta_abs_2 = 0;

	for (;theta_dof < theta_dof_end; ++gdof
							 	   , ++subdof
							 	   , ++theta_dof)
	{
		// if not released dof
		if (*subdof == 0) {
			// if constrained main dof
			if (*gdof == 0) {
				*theta_dof = 0.0;
				continue;
			}
			// main dof
			*theta_dof = q[*gdof-1];
		} else {
			// released dof
			*theta_dof = q[*subdof-1];
			q[*subdof-1] = 0;
		}
		theta_abs_2 += (*theta_dof)*(*theta_dof);
	}

	if (theta_abs_2 > 0) {
		// begin of the sub node Rsum
		// std::cout << "here" << std::endl;
		math::Slice Rsumi_sub(Rsumi,Rsumi+3);
		// Rsumi_sub = math::dot(math::rotation_tensor(theta), Rsumi_sub);
		math::rotation_tensor(theta,temp_rotTensor);
		math::fill(temp_Rsumi.begin(),temp_Rsumi.end(),0.0);
		math::dot(temp_rotTensor,Rsumi_sub,temp_Rsumi);
		math::swap(Rsumi_sub,temp_Rsumi);
	}
}

void AnalysisTraits::update_Rsumi_releases(math::vector<double>& q
								, typename math::vector_t<double,2>::iterator Rsumi
								, typename math::vector<size_t>::const_iterator gdof
								, 		   math::vector<double>& theta
								,		   math::vector_t<double,2>& temp_rotTensor
								,		   math::vector_t<double,2>& temp_Rsumi)
{
	auto theta_dof = theta.begin()
		,theta_dof_end = theta.end();
	double theta_abs_2 = 0;

	for (;theta_dof < theta_dof_end; ++gdof
							 	   , ++theta_dof)
	{
		if (*gdof == 0) {
			*theta_dof = 0.0;
			continue;
		}
		
		*theta_dof = q[*gdof-1]; // 
		q[*gdof-1] = 0;          //
		theta_abs_2 += (*theta_dof)*(*theta_dof);
	}

	if (theta_abs_2 > 0) {
		// begin of the sub node Rsum
		// std::cout << "here" << std::endl;
		math::Slice Rsumi_sub(Rsumi,Rsumi+3);
		// Rsumi_sub = math::dot(math::rotation_tensor(theta), Rsumi_sub);
		math::rotation_tensor(theta,temp_rotTensor);
		math::fill(temp_Rsumi.begin(),temp_Rsumi.end(),0.0);
		math::dot(temp_rotTensor,Rsumi_sub,temp_Rsumi);
		math::swap(Rsumi_sub,temp_Rsumi);
	}
}

void AnalysisTraits::update_Rsumi(math::vector<double>& q
								, typename math::vector_t<double,3>::iterator Rsum_node
								, typename math::vector<size_t>::const_iterator gdof
								, 		   math::vector<double>& theta
								,		   math::vector_t<double,2>& temp_rotTensor
								,		   math::vector_t<double,2>& temp_Rsumi)
{
	auto theta_dof = theta.begin()
		,theta_dof_end = theta.end();
	double theta_abs_2 = 0;

	for (;theta_dof < theta_dof_end; ++gdof
							 	   , ++theta_dof)
	{
		if (*gdof == 0) {
			*theta_dof = 0.0;
			continue;
		}
		
		*theta_dof = q[*gdof-1]; // 
		q[*gdof-1] = 0;          //
		theta_abs_2 += (*theta_dof)*(*theta_dof);
	}

	if (theta_abs_2 > 0) {
		// begin of the sub node Rsum
		// std::cout << "here" << std::endl;
		// Rsumi_sub = math::dot(math::rotation_tensor(theta), Rsumi_sub);
		math::rotation_tensor(theta,temp_rotTensor);
		math::fill(temp_Rsumi.begin(),temp_Rsumi.end(),0.0);
		math::dot(temp_rotTensor,*Rsum_node,temp_Rsumi);
		math::swap(*Rsum_node,temp_Rsumi);
	}
}

math::vector_t<double,3> AnalysisTraits::setup_Rsum(const Model& model) {
	math::vector_t<double,3> Rsum(model.nodes_info.size(),math::eye<double>(3));

#ifdef BEAM_RELEASES
	size_t num_subnodes; // number of released elements in node
	size_t subnode;
	for (const NodeAllReleases& node : model.nodes_releases) {
		num_subnodes = node.dofs.size()/BaseNode::DOFS;
		Rsum[node.node_globalid].resize(3*(num_subnodes+1));
		for (subnode = 0; subnode < num_subnodes; ++subnode) {
			Rsum[node.node_globalid][(subnode+1)*3  ] = {1.0, 0.0, 0.0};
			Rsum[node.node_globalid][(subnode+1)*3+1] = {0.0, 1.0, 0.0};
			Rsum[node.node_globalid][(subnode+1)*3+2] = {0.0, 0.0, 1.0};
		}
	}
#endif

	return Rsum;
}




} // namespace fem