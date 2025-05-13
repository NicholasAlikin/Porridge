#pragma once

// std headers

// local headers
#include "linalg.hpp"

#include "basic_system.hpp"

#include "fem.hpp"
#include "analyzes.hpp"

// other headers


namespace fem::npath {
// namespace npath {
/*Main euqation:
    h(q,l) = l*g - f(q) = 0
where
    q - displacement vector
    l - load factor ()
    g - unitial external load vector
*/
class StaticNL: public ::npath::BasicSystem<> {
private:
public:
    
    Model model;
    Assemble assemble;
    
    math::vector<double> loadExt;        // external load, |loadExt| == 1
    math::vector<double> loadInt;
    math::vector_t<double,3> Rsum;
    math::vector<size_t> ynorm_pos; // displacement vector components, which need to print

    double loadExt_norm;
    
    math::vector_t<double,4> Rsum_continuation; // Rsum values for Continuation::continuation.process() for each curve point

    math::vector<double> loadExt_unit;
private:
    math::vector<double> temp_Dy;
    math::vector<double> temp_theta;
	math::vector_t<double,2> temp_rotTensor;
	math::vector_t<double,2> temp_Rsumi;

    bool incremental_large_rotation;
public:

    StaticNL() = default;
    
    StaticNL(const StaticNL& other);
    
    StaticNL(StaticNL&& other);
    
    StaticNL(const Model& model
            ,const Assemble& assemble
            ,const math::vector<size_t>& ynorm_pos
            ,bool incremental_large_rotation = true);
    
    void initialization();
    
    // void process_initialization();
    void continuation_initialization(math::vector_t<double,2>& jac
                                    ,math::vector<double>& ynorm);


    void process_total_increment(math::vector<double>& y
                               , math::vector<double>& dy);

    void process_total_increment(math::vector<double>& Dy
                               , math::vector<double>& dy
                         , const math::vector<double>& predictor
                        ,                     double   ds);


    void system_response(math::vector<double>&     fun
                       , math::vector_t<double,2>& jac
                       , math::vector<double>&     y);
    void system_response_extended(math::vector<double>& fun
                                , math::vector_t<double,2>& jac
                                , math::vector<double>& y);
    

    void calculate_response_norm(const math::vector<double>& y
                                     , math::vector<double>& ynorm);
    size_t response_norm_size();
    void solution_not_found();
    void step_back(double new_step);

    double fun_norm(const math::vector<double>& fun, const math::vector<double>& y);



    void linear_like_solve(math::vector_t<double,2>& stif, math::vector<double>& y);
    void update_Rsum(math::vector<double>& q);
    void update_loadExt_vector(double new_loadExt_norm);

private:
    void prepare_fun_jac(math::vector<double>& fun
                        ,math::vector_t<double,2>& jac);
    void zeros_stiffness_matrix(math::vector_t<double,2>& stif);

    void do_assemble(math::vector_t<double,2>& stif, const math::vector<double>& q);
};



} // namespace fem::npath