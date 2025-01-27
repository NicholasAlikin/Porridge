#pragma once

// std headers

// local headers
#include "linalg.hpp"

#include "basic_system.hpp"

#include "fem.hpp"
#include "analyzes.hpp"

// other headers


namespace fem::npath {

#if 0
class HBM: public ::npath::BasicSystem {


private:
    math::vector<double> func_ext;
    math::vector<double> func_inerc;
    math::vector<double> func_damp;
    math::vector<double> func_inner;

};
#endif


class ShootingMethod: public ::npath::BasicSystem {
public:
    Model model;
    Assemble assemble;

    math::vector<double> load_ext;       // external load, |loadExt| == 1
    math::vector<double> load_inner;     // inner forces (internal)
    math::vector<double> load_inert;     // inertia forces
    math::vector_t<double,3> Rsum;
    math::vector<size_t> ynorm_pos; // displacement vector components, which need to print

    math::vector_t<double,4> Rsum_continuation; // Rsum values for Continuation::continuation.process() for each curve point

public:
    /* necessary BaiscSystem methods */
    void system_response(math::vector<double>&     fun
                       , math::vector_t<double,2>& jac
                       , math::vector<double>&     y);
    void system_response_extended(math::vector<double>&     fun
                                , math::vector_t<double,2>& jac
                                , math::vector<double>&     y);


    void calculate_response_norm(const math::vector<double>& y
                                     , math::vector<double>& ynorm);
    size_t response_norm_size();
    
    /* not necessary BaiscSystem methods */
    void solution_not_found();
    void step_back(double new_step);

};

} // namespace fem::npath