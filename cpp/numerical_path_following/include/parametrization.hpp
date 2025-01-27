#pragma once

// std headers

// local headers
#include "linalg.hpp"
#include "basic_system.hpp"

#include "numerical_continuation_base.hpp"

// other headers


/*Methods of curve parametrization*/

namespace npath {

struct BasePathFollower {
    static const bool do_correction = true;
    static const bool do_continuation = false;
    static const bool need_extended_system_response = false;
protected:
    size_t sys_size;
public:
    BasePathFollower(size_t system_size);
    ~BasePathFollower() = default;

    size_t system_size() const;
    size_t fun_size() const;
    std::array<size_t,2> jac_size() const;
    

    void process_iteration(math::vector<double>&     fun
                         , math::vector_t<double,2>& jac
                         , math::vector<double>&     dy);

    void process_total_increment(math::vector<double>& y
                         , const math::vector<double>& dy);
    void process_total_increment(math::vector<double>& Dy
                         , const math::vector<double>& dy
                         , const math::vector<double>& predictor
                         ,                    double   ds) = delete;
    
    void system_response_extended(math::vector<double>&     fun
                          ,       math::vector_t<double,2>& jac
                          , const math::vector<double>&     Dy
                          ,                    double       ds) = delete;
};

struct NormalFlow: BasePathFollower {
    static const bool do_correction = true;
    static const bool do_continuation = true;
    static const bool need_extended_system_response = true;
protected:
    math::vector_t<double,2>  temp_AAT;
    math::vector_t<double,2>  temp_L; // for internal calculations, to not allocate memory for each iteration
    math::vector_t<double,1>  temp_D;

public:
    NormalFlow(size_t system_size);
    
    
    std::array<size_t,2> jac_size() const;

    void process_iteration(math::vector<double>&     fun
                         , math::vector_t<double,2>& jac
                         , math::vector<double>&     dy);

    void process_total_increment(math::vector<double>& y
                         , const math::vector<double>& dy);

    void process_total_increment(math::vector<double>& Dy
                         , const math::vector<double>& dy
                         , const math::vector<double>& predictor
                         ,                    double   ds);
    
    void system_response_extended(math::vector<double>&     fun
                          ,       math::vector_t<double,2>& jac     
                          , const math::vector<double>&     Dy
                          ,                    double       ds);
};


struct NormalFlowNoStepCorrection: NormalFlow {
    using NormalFlow::NormalFlow;
    using NormalFlow::process_total_increment;

    void process_total_increment(math::vector<double>& Dy
                         , const math::vector<double>& dy
                         , const math::vector<double>& predictor
                         ,                    double   ds);
};

struct NormalFlowNoStepCorrectionStructural: NormalFlowNoStepCorrection {
    using NormalFlowNoStepCorrection::NormalFlowNoStepCorrection;
    void process_iteration(math::vector<double>&     fun
                         , math::vector_t<double,2>& jac
                         , math::vector<double>&     dy);

protected:
    math::vector_t<double,2> invW = math::eye<double>(sys_size);
public:
    void set_inverse_weight_matrix(const math::vector_t<double,2>& invW);
};

struct ArcLength: BasePathFollower {
    static const bool do_correction = false;
    static const bool do_continuation = true;
    static const bool need_extended_system_response = true;

    size_t fun_size() const;
    std::array<size_t,2> jac_size() const;

    void process_total_increment(math::vector<double>& y
                         , const math::vector<double>& dy) = delete;
    void process_total_increment(math::vector<double>& Dy
                         , const math::vector<double>& dy
                         , const math::vector<double>& predictor
                         ,                    double   ds);
    
    void system_response_extended(math::vector<double>&     fun
                          ,       math::vector_t<double,2>& jac
                          , const math::vector<double>&     Dy
                          ,                    double       ds);
};

} // namespace npath