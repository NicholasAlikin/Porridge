#include "staticNL.hpp"

namespace fem::npath {

StaticNL::StaticNL(const StaticNL& other)
        : model(other.model)
        , assemble(other.assemble)
        
        , loadExt(other.loadExt)
        , loadInt(other.loadInt)
        , Rsum(other.Rsum)
        , ynorm_pos(other.ynorm_pos)
        
        , loadExt_norm(other.loadExt_norm)
        
        , Rsum_continuation(other.Rsum_continuation)
        
        , temp_Dy(other.temp_Dy)
        , temp_theta(other.temp_theta)
        , temp_rotTensor(other.temp_rotTensor)
        , temp_Rsumi(other.temp_Rsumi)

        ,incremental_large_rotation(other.incremental_large_rotation)
{}

StaticNL::StaticNL(StaticNL&& other)
        : model(std::move(other.model))
        , assemble(std::move(other.assemble))
        
        , loadExt(std::move(other.loadExt))
        , loadInt(other.loadInt)
        , Rsum(std::move(other.Rsum))
        , ynorm_pos(std::move(other.ynorm_pos))
        
        , loadExt_norm(other.loadExt_norm)
        
        , Rsum_continuation(std::move(other.Rsum_continuation))
        
        , temp_Dy(std::move(other.temp_Dy))
        , temp_theta(std::move(other.temp_theta))
        , temp_rotTensor(std::move(other.temp_rotTensor))
        , temp_Rsumi(std::move(other.temp_Rsumi))

        ,incremental_large_rotation(other.incremental_large_rotation)
{}

StaticNL::StaticNL(const Model& model
                  ,const Assemble& assemble
                  ,const math::vector<size_t>& ynorm_pos
                  ,bool incremental_large_rotation)
        : model(model)
        , assemble(assemble)
        , ynorm_pos(ynorm_pos)
        , incremental_large_rotation(incremental_large_rotation)
{
    initialization();
}

void StaticNL::initialization() {
    AnalysisTraits::setup_analysis(model,assemble);

    loadExt = math::zeros<double>(assemble.ndofs);
    ModelTraits::store_load_vector(model,assemble,loadExt);
    loadExt_norm = norm(loadExt);
    
    loadInt = math::zeros<double>(assemble.ndofs);
    
    
    if (incremental_large_rotation) {
        temp_Dy = math::zeros<double>(assemble.ndofs+1);
        
        Rsum = AnalysisTraits::setup_Rsum(model);
        temp_Rsumi = math::zeros<double>(3,3);
        temp_theta = math::zeros<double>(3);
        temp_rotTensor = math::zeros<double>(3,3);
    }
    
    loadExt_unit = loadExt/loadExt_norm;
}



void StaticNL::continuation_initialization(math::vector_t<double,2>& jac
                                          ,math::vector<double>& ynorm) {
    Rsum_continuation = math::vector_t<double,4>();
    
    // last jac column is constant and equal to negative unary external load
    size_t last = assemble.ndofs;
    auto jac_row = jac.begin();
    auto load_ext = loadExt.begin(), load_ext_end = loadExt.end();
    while (load_ext < load_ext_end) {
        (*jac_row)[last] = -(*load_ext)/loadExt_norm;
        ++jac_row; ++load_ext;
    }

    // set ynorm size
    ynorm.resize(ynorm_pos.size());
}


void StaticNL::process_total_increment(math::vector<double>& y
                                     , math::vector<double>& dy) {
    update_Rsum(dy);
}

void StaticNL::process_total_increment(math::vector<double>& Dy
                                     , math::vector<double>& dy
                               , const math::vector<double>& predictor
                               ,                    double   ds) {
    if (!incremental_large_rotation) return;

#if 1
    double dot_pre_Dy = math::dot(predictor,Dy)
         , dot_pre_dy = math::dot(predictor,dy);
    
    if (dot_pre_Dy+dot_pre_dy <= 0) {
        double backward_coef = -dot_pre_Dy/dot_pre_dy * ::npath::PathFollowing::track_backward;
        Dy + dy*backward_coef;
    } else {
        math::sum(Dy.begin(),dy.begin(),temp_Dy.begin(),temp_Dy.end());
    }
#else
    temp_Dy = Dy + dy;
#endif

    dy = temp_Dy*(ds/math::norm(temp_Dy)) - Dy;

    update_Rsum(dy);
}


void StaticNL::system_response(math::vector<double>& fun
                             , math::vector_t<double,2>& jac
                             , math::vector<double>& y)
{   
    update_loadExt_vector(y.last());
    
    // assembleNL uses y by only indeses of q (displacement vector)
    prepare_fun_jac(fun,jac);
    do_assemble(jac,y);
    // Calculate fun = loadExt - loadInt
    math::sub(loadExt.begin(), loadInt.begin(), fun.begin(), fun.begin()+loadExt.size());

    // std::cout << "# stif[0] = " << math::Slice(jac[0].begin(),jac[0].end()) << "\n";
    
}

void StaticNL::system_response_extended(math::vector<double>& fun
                                      , math::vector_t<double,2>& jac
                                      , math::vector<double>& y)
{   
    system_response(fun,jac,y);
}


void StaticNL::calculate_response_norm(const math::vector<double>& y
                                           , math::vector<double>& ynorm)
{
    auto it_ynorm = ynorm.begin()
        ,it_norm_end = ynorm.end();
    auto it_ynorm_pos = ynorm_pos.begin();

    while (it_ynorm < it_norm_end) {
        *it_ynorm = y[*it_ynorm_pos];
        ++it_ynorm; ++it_ynorm_pos;
    }

    Rsum_continuation.push_back(Rsum);
}

size_t StaticNL::response_norm_size() {
    return ynorm_pos.size();
}




void StaticNL::linear_like_solve(math::vector_t<double,2>& stif, math::vector<double>& y) {
    update_loadExt_vector(y.last());
    
    do_assemble(stif,y);
    
    math::solve_ldlt(stif,loadExt,y);
}

void StaticNL::update_Rsum(math::vector<double>& q) {
    if (!incremental_large_rotation) return;

    AnalysisTraits::update_Rsum(model,assemble,Rsum,q
                            ,temp_theta,temp_rotTensor,temp_Rsumi);
}

void StaticNL::update_loadExt_vector(double new_loadExt_norm) {
	auto load_node     = model.loads_info.begin()
		,load_node_end = model.loads_info.end();
		
	decltype(model.loads_info[0].dofs.begin()) load_dof, load_dof_end;
	
	decltype(assemble.nodes_dofs.begin()) gnode;
	decltype(assemble.nodes_dofs[0].begin()) gdof;
	double loadExt_norm_frac = new_loadExt_norm/loadExt_norm;
	while (load_node != load_node_end) {
		load_dof = load_node->dofs.begin();
		load_dof_end = load_node->dofs.end();
		gdof = assemble.nodes_dofs[load_node->id].begin();
		while (load_dof != load_dof_end) {
			if (*gdof > 0) {
                *load_dof *= loadExt_norm_frac;
				loadExt[*gdof-1] = *load_dof;
			}
			++load_dof; ++gdof;
		}
		++load_node;
	}
    loadExt_norm = new_loadExt_norm;

}

void StaticNL::solution_not_found() {
    if (incremental_large_rotation)
        std::copy(Rsum_continuation.last().begin(),Rsum_continuation.last().end(),Rsum.begin());
    // Rsum = Rsum_continuation.last();
}

void StaticNL::step_back(double) {
    if (incremental_large_rotation) {
        Rsum_continuation.erase(Rsum_continuation.end()-1);
        std::copy(Rsum_continuation.last().begin(),Rsum_continuation.last().end(),Rsum.begin());
    }
    // Rsum = *(Rsum_continuation.end()-1);
    // std::cout << "# SRsum = \n" << Rsum << std::endl;
    // throw 1;
}

double StaticNL::fun_norm(const math::vector<double>& fun, const math::vector<double>& y) {
    return math::norm(fun) / y.last();
}


void StaticNL::prepare_fun_jac(math::vector<double>& fun
                              ,math::vector_t<double,2>& jac) {
    math::fill(loadInt.begin(),loadInt.end(), 0.0);
    
    zeros_stiffness_matrix(jac);
}

void StaticNL::zeros_stiffness_matrix(math::vector_t<double,2>& stif) {
    for (auto& row: stif) {
        math::fill(row.begin(),row.begin()+assemble.ndofs, 0.0);
    }
    // stif = zeros<double>(stif);

    // // last jac column is constant and equal to negative unary external load
    // size_t last = assemble.ndofs;
    // auto jac_row = stif.begin();
    // auto load_ext = loadExt_unit.begin(), load_ext_end = loadExt_unit.end();
    // while (load_ext < load_ext_end) {
    //     (*jac_row)[last] = -(*load_ext);
    //     ++jac_row; ++load_ext;
    // }
}


void StaticNL::do_assemble(math::vector_t<double,2>& stif, const math::vector<double>& q) {
    if (incremental_large_rotation) {
        ModelTraits::assemble(model,assemble,stif,loadInt,q,Rsum,&BaseElement::tangentStiffness_innerLoad);
    }
    else {
        ModelTraits::assemble(model,assemble,stif,loadInt,q,&BaseElement::tangentStiffness_innerLoad);
    }
}

} // namespace fem::npath