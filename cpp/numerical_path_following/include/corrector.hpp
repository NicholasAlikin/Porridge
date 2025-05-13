/*
Corrector algorithms:
    o default: Newton-Raphson algorithm
    o normal-flow: Gauss-Newton algorithm
                   using Moore-Penrose inverse
                   (pseudoinverse) matrix
    o Psevdo-arc-lenght algorithms:
        - sphere scheme,
        - orthogonal scheme.
*/
#pragma once

// std headers

// local headers
#include "numerical_continuation_base.hpp"
#include "linalg.hpp"

// other headers


namespace npath { 

/*
BasicSystem is neaded only to calculate system response (fun and jac).
        It rewrites values to fun and jac, which correspond `y` value.

PathFollower is neaded to do sub iteration (solve linear equation),
        it also store temporal VectorLike objects which are neaded while solving linear equations
        It rewrites linear equation answer `dy`.

Corrector is merge BasicSystem and PathFollower. It calls all methods.
        It store temporal VectorLike objects which are neaded to do iteration:
        Every process mast starts with zeros neaded temporal VectorLike objects (`Dy`).

*/
template <typename BasicSystem, typename PathFollower>
class Corrector {
public:
    
    BasicSystem basic_system;
    PathFollower path_follower;
    using jac_t = PathFollower::jac_t;
// private:
    jac_t jac;
    math::vector<double>     fun;
    math::vector<double>     dy;
    math::vector<double>     Dy;
public:
    size_t iter = 0, global_iter = 0;
    int exitflag = 0;
    double _dy_norm = 0, _fun_norm = 0;
    double _epsx = 1e-5, _epsf = 1e-5;
    double _dx = 1e-5;
    size_t _max_iter = 8;

public:

    Corrector() = default;
    Corrector(const BasicSystem& basic_system, const PathFollower& path_follower);

    void process(math::vector<double>& y);

    void process(math::vector<double>& y
         , const math::vector<double>& predictor
         , const math::vector<double>& previous
                            , double   ds);

    void process_iteration(math::vector<double>& y);
    void process_iteration(math::vector<double>& y
                   , const math::vector<double>& predictor
                   , const math::vector<double>& previous
                                      , double ds);

    void process_initialization();
    void continuation_initialization(math::vector<double>& ynorm);
    void process_exitflag(const math::vector<double>& dy
                        , const math::vector<double>& fun
                        , const math::vector<double>& y);
    void printiter() const;

    void calculate_response_norm(const math::vector<double>& y
                                     , math::vector<double>& ynorm);
    size_t response_norm_size();
    void solution_not_found();
    void step_back(double new_step);
    
};

template <typename BasicSystem, typename PathFollower>
Corrector<BasicSystem,PathFollower>::Corrector(const BasicSystem& basic_system
                                             , const PathFollower& path_follower)
        : basic_system(basic_system)
        , path_follower(path_follower)
        , jac(path_follower.jac_init())
        , fun(math::zeros<double>(path_follower.fun_size()))
        , dy (math::zeros<double>(path_follower.system_size()))
        , Dy (math::zeros<double>(path_follower.system_size()))
{}


template <typename BasicSystem, typename PathFollower>
Corrector(const BasicSystem&, const PathFollower&)
        -> Corrector<BasicSystem,PathFollower>;


template <typename BasicSystem, typename PathFollower>
void Corrector<BasicSystem,PathFollower>::process(math::vector<double>& y) {
    process_initialization();
    math::fill(Dy.begin(),Dy.end(),0.0);

    process_iteration(y);
    while (exitflag == EXITFLAG::OK) {
        process_iteration(y);
    }
}

template <typename BasicSystem, typename PathFollower>
void Corrector<BasicSystem,PathFollower>::process(math::vector<double>& y
                                           , const math::vector<double>& predictor
                                           , const math::vector<double>& previous
                                           ,                    double   ds) {
    /* Its assumed, that `y` was calculated by this method,
    so it is equal to previos */
    process_initialization();
    
    y += predictor;
    std::copy(predictor.begin(),predictor.end(),Dy.begin());

    process_iteration(y,predictor,previous,ds);
    while (exitflag == EXITFLAG::OK) {
        process_iteration(y,predictor,previous,ds);
    }
    global_iter += iter;
}

template <typename BasicSystem, typename PathFollower>
void Corrector<BasicSystem,PathFollower>::process_initialization() {
    iter = 0;
}

template <typename BasicSystem, typename PathFollower>
void Corrector<BasicSystem,PathFollower>::continuation_initialization(math::vector<double>& ynorm) {
    global_iter = 0;
    
    basic_system.continuation_initialization(jac,ynorm);
}

template <typename BasicSystem, typename PathFollower>
void Corrector<BasicSystem,PathFollower>::process_iteration(math::vector<double>& y) {
    // 0. increment iteration counter
    ++iter;
    // 1. calculate system response `fun` and `jac` for point `y`
    if constexpr (PathFollower::need_extended_system_response) {
        basic_system.system_response_extended(fun,jac,y);
    } else {
        basic_system.system_response(fun,jac,y);
    }
    // 2. calculate increment dy
    path_follower.process_iteration(fun,jac,dy);

    // 3. calcualte exitflag
    process_exitflag(dy,fun,y);

    // 4. calculate total increment for basic_system state ...
    basic_system.process_total_increment(y,dy);
    // ... and using parametrization method
    path_follower.process_total_increment(y,dy);

    // 5. print iteration info
    printiter();
}

template <typename BasicSystem, typename PathFollower>
void Corrector<BasicSystem,PathFollower>::process_iteration(math::vector<double>& y
                                                    , const math::vector<double>& predictor
                                                    , const math::vector<double>& previous
                                                                       , double ds)
{
    // 0. increment iteration counter
    ++iter;
    // 1. calculate system response -fun and jac for point y
    if constexpr (PathFollower::need_extended_system_response) {
        basic_system.system_response_extended(fun,jac,y);
    } else {
        basic_system.system_response(fun,jac,y);
    }
    // 1.1. sub-calculate system response using parametrization method
    path_follower.system_response_extended(fun,jac,Dy,ds);

    // 2. calculate increment dy
    path_follower.process_iteration(fun,jac,dy);

    // 3. calcualte exitflag
    process_exitflag(dy,fun,y);

    // 4. calculate total increment `Dy` for basic_system state ...
    basic_system.process_total_increment(Dy,dy,predictor,ds);
    // 4.1. and using parametrization method
    path_follower.process_total_increment(Dy,dy,predictor,ds);
    
    // 5. calculate new point `y = Dy + predictor`
    math::sum(Dy.begin(),previous.begin(),y.begin(),y.end());
    
    // 6. print iteration info
    printiter();
}


template <typename BasicSystem, typename PathFollower>
void Corrector<BasicSystem,PathFollower>::process_exitflag(const math::vector<double>& dy
                                                         , const math::vector<double>& fun
                                                         , const math::vector<double>& y) {
    _dy_norm = math::norm(dy);
    _fun_norm = basic_system.fun_norm(fun,y);

    if (_dy_norm < _epsx && _fun_norm < _epsf) {
        exitflag = EXITFLAG::NORM_VAR_AND_FUN;
    } else if (iter >= _max_iter) {
        exitflag = EXITFLAG::MAX_ITER;
    } else {
        exitflag = EXITFLAG::OK;
    }
}

template <typename BasicSystem, typename PathFollower>
void Corrector<BasicSystem,PathFollower>::printiter() const {
#ifdef CORRECTOR_PRINTITER
    std::cout << "#\t Corrector Iter: " << iter
              << ", |f|: " << _fun_norm
              << ", |dy|: " << _dy_norm
              << '\n';
#endif
}

template <typename BasicSystem, typename PathFollower>
void Corrector<BasicSystem,PathFollower>::calculate_response_norm(const math::vector<double>& y
                                                                      , math::vector<double>& ynorm) {
    basic_system.calculate_response_norm(y,ynorm);
}
template <typename BasicSystem, typename PathFollower>
size_t Corrector<BasicSystem,PathFollower>::response_norm_size() {
    return basic_system.response_norm_size();
}

template <typename BasicSystem, typename PathFollower>
void Corrector<BasicSystem,PathFollower>::solution_not_found() {
    basic_system.solution_not_found();
}

template <typename BasicSystem, typename PathFollower>
void Corrector<BasicSystem,PathFollower>::step_back(double new_step) {
    basic_system.step_back(new_step);
}




} // namespace npath