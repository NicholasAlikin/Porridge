#include "parametrization.hpp"

namespace npath {

/*------------
BasePathFollower
------------*/
BasePathFollower::BasePathFollower(size_t system_size)
        : sys_size(system_size)
{}

size_t BasePathFollower::system_size() const {
    return sys_size;
}

size_t BasePathFollower::fun_size() const {
    return sys_size-1;
}
std::array<size_t,2> BasePathFollower::jac_size() const {
    return {sys_size-1,sys_size-1};
}

void BasePathFollower::process_iteration(math::vector<double>&     fun
                                       , math::vector_t<double,2>& jac
                                       , math::vector<double>&     dy) {
    math::solve_lu(jac,fun,dy,sys_size-1);
}

void BasePathFollower::process_total_increment(math::vector<double>& y
                                       , const math::vector<double>& dy) {
    math::Slice sly(y.begin(), y.end()-1);
    sly += dy;
}


/*------------
NormalFlow
------------*/
NormalFlow::NormalFlow(size_t system_size)
        : BasePathFollower(system_size)
        , temp_AAT(math::zeros<double>(system_size-1,system_size-1))
        , temp_L(math::zeros<double>(system_size-1,system_size-1))
        , temp_D(math::zeros<double>(system_size-1))
{}

std::array<size_t,2> NormalFlow::jac_size() const {
    return {sys_size-1,sys_size};
}

void NormalFlow::process_iteration(math::vector<double>&     fun
                                 , math::vector_t<double,2>& jac
                                 , math::vector<double>&     dy) {
    math::psolve(jac,fun,dy,temp_AAT,temp_L,temp_D);
}

void NormalFlow::process_total_increment(math::vector<double>& y
                                 , const math::vector<double>& dy) {
    y += dy;
}

void NormalFlow::process_total_increment(math::vector<double>& Dy
                                 , const math::vector<double>& dy
                                 , const math::vector<double>& predictor
                                                    , double   ds) {
    // backward moment
    double dot_pre_Dy = math::dot(predictor,Dy)
         , dot_pre_dy = math::dot(predictor,dy);
    if (dot_pre_Dy+dot_pre_dy <= 0) {
        double backward_coef = -dot_pre_Dy/dot_pre_dy * PathFollowing::track_backward;
        Dy += dy*backward_coef;
    } else {
        Dy += dy;
    }

    // step length
    Dy *= ds/math::norm(Dy);
}


void NormalFlow::system_response_extended(math::vector<double>&
                                        , math::vector_t<double,2>&
                                  , const math::vector<double>&
                                                    , double ds)
{}



/*------------
NormalFlowNoStepCorrection
------------*/
void NormalFlowNoStepCorrection::process_total_increment(math::vector<double>& Dy
                                                 , const math::vector<double>& dy
                                                 , const math::vector<double>&
                                                                    , double) {
    Dy += dy;
}

/*------------
NormalFlowNoStepCorrectionStructural
------------*/

void NormalFlowNoStepCorrectionStructural::process_iteration(math::vector<double>&     fun
                                 , math::vector_t<double,2>& jac
                                 , math::vector<double>&     dy) {
    math::psolve_weight(jac,fun,invW,dy);
}

void NormalFlowNoStepCorrectionStructural::set_inverse_weight_matrix(const math::vector_t<double, 2>& invW) {
    this->invW = invW;
}

/*------------
ArcLength
------------*/

size_t ArcLength::fun_size() const {
    return sys_size;
}
std::array<size_t,2> ArcLength::jac_size() const {
    return {sys_size,sys_size};
}

void ArcLength::process_total_increment(math::vector<double>& Dy
                                , const math::vector<double>& dy
                                , const math::vector<double>& predictor
                                                   , double) {
    // backward moment
    double dot_pre_Dy = math::dot(predictor,Dy)
         , dot_pre_dy = math::dot(predictor,dy);
    if (dot_pre_Dy+dot_pre_dy <= 0) {
        double backward_coef = -dot_pre_Dy/dot_pre_dy * PathFollowing::track_backward;
        Dy += dy*backward_coef;
    } else {
        Dy += dy;
    }
}

void ArcLength::system_response_extended(math::vector<double>&     fun
                                 ,       math::vector_t<double,2>& jac
                                 , const math::vector<double>&     Dy
                                 ,                    double       ds) {
    fun.last() = -(dot(Dy,Dy) - ds*ds);
    jac.last() = 2*Dy;
}


} // namespace npath