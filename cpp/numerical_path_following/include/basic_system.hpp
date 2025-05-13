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
template <typename jac_t = math::vector_t<double,2>>
class BasicSystem {
public:
    
    BasicSystem() = default;
    ~BasicSystem() = default;

public:

    void continuation_initialization(jac_t& jac
                                    ,math::vector<double>& ynorm);  


    void process_total_increment(math::vector<double>& y
                               , math::vector<double>& dy);

    void process_total_increment(math::vector<double>& Dy
                               , math::vector<double>& dy
                         , const math::vector<double>& predictor
                         ,                    double   ds);


    void system_response(math::vector<double>&     fun
                       , jac_t& jac
                       , math::vector<double>&     y) = delete;

    void system_response_extended(math::vector<double>&     fun
                                , jac_t& jac
                                , math::vector<double>&     y) = delete;


    void calculate_response_norm(const math::vector<double>& y
                                     , math::vector<double>& ynorm) = delete;
    size_t response_norm_size() = delete;
    void solution_not_found();
    void step_back(double new_step);
    double fun_norm(const math::vector<double>& fun, const math::vector<double>& y);
};

template <typename jac_t>
void BasicSystem<jac_t>::continuation_initialization(jac_t& jac
                                             ,math::vector<double>& ynorm) {
      // ynorm.resize(response_norm_size());
}


template <typename jac_t>
void BasicSystem<jac_t>::process_total_increment(math::vector<double>&
                                        , math::vector<double>&) {
}


template <typename jac_t>
void BasicSystem<jac_t>::process_total_increment(math::vector<double>&
                                        , math::vector<double>&
                                  , const math::vector<double>&
                                  ,                    double) {
}


template <typename jac_t>
void BasicSystem<jac_t>::solution_not_found() {
}


template <typename jac_t>
void BasicSystem<jac_t>::step_back(double) {
}


template <typename jac_t>
double BasicSystem<jac_t>::fun_norm(const math::vector<double>& fun,
                                    const math::vector<double>&) {
    return math::norm(fun);
}

} // namespace npath