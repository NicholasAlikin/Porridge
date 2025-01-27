#pragma once

// std headers

// local headers
#include "linalg.hpp"
#include "numerical_continuation_base.hpp"
#include "corrector.hpp"
// other headers


namespace npath {

/*System to define underdetermined system of equation,
which implicitly define curve.
For example
    - for dynamic analysis: harmonic balance method, shooting method;
    - for static analysis : quilibrium equation. */
class BasicSystem {
public:
    
    BasicSystem() = default;
    ~BasicSystem() = default;

public:

    void continuation_initialization(math::vector_t<double,2>& jac
                                    ,math::vector<double>& ynorm);    


    void process_total_increment(math::vector<double>& y
                               , math::vector<double>& dy);

    void process_total_increment(math::vector<double>& Dy
                               , math::vector<double>& dy
                         , const math::vector<double>& predictor
                         ,                    double   ds);


    void system_response(math::vector<double>&     fun
                       , math::vector_t<double,2>& jac
                       , math::vector<double>&     y) = delete;

    void system_response_extended(math::vector<double>&     fun
                                , math::vector_t<double,2>& jac
                                , math::vector<double>&     y) = delete;


    void calculate_response_norm(const math::vector<double>& y
                                     , math::vector<double>& ynorm) = delete;
    size_t response_norm_size() = delete;
    void solution_not_found();
    void step_back(double new_step);
};

} // namespace npath