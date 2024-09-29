#pragma once

#include "numerical_continuation_base.hpp"
#include "linalg.hpp"
#include "dft.hpp"
#include "predictor.hpp"
#include "corrector.hpp"


namespace math {

template <typename Method, typename Pred>
class Continuation {

public:
    
    Method method;
    Pred predictor;
    double _ds;
    int exitflag;
    size_t processiter;
    size_t _successful_steps;
    size_t _count_not_correct_solution;



    Continuation() = default;
    Continuation(const Method& method, const Pred& predictor);

    template <typename ParMethod>
    requires ParMethod::do_continuation
    void process(const Vector_t& x0, double param_start, double param_end, double ds0
        , bool if_printiter, double arclen_min=1e-4, double arclen_max=1.0, double arclen_inc=2.0
        , size_t successful_steps_max=10);
    
    void process_sub_iteration(const Vector_t& y, Matrix_t& Y, Matrix_t& Ynorm, double param_end, bool if_printiter);
    
    void process_initialization(double param_start, double param_end, double ds0);

    template <typename ParMethod>
    Vector_t zeros_iteration(const Vector_t& x0,double param_start);
    
    void process_exitflag(double param, double param_end);
    void printiter(const Vector_t& ynorm, double param, bool if_printiter);
    Vector_t calculate_response_norm(const Vector_t& y, size_t ndof);
    void process_end_message();
    void step_update(double arclen_min, double arclen_max, double arclen_inc
                    ,size_t successful_steps_max);
    bool decrease_step(double arclen_min, double arclen_inc);
    bool increase_step(double arclen_max, double arclen_inc, size_t successful_steps_max);
    bool is_correct_solution(const Vector_t& y, const Matrix_t& Y, double arclen_min, double arclen_inc);
};


template <typename Method, typename Pred>
Continuation<Method,Pred>::Continuation(const Method& method, const Pred& predictor)
        : method(method), predictor(predictor)
{

}



template <typename Method, typename Pred>
template <typename ParMethod>
requires ParMethod::do_continuation
void Continuation<Method,Pred>::process(const Vector_t& x0, double param_start, double param_end, double ds0
                , bool if_printiter, double arclen_min, double arclen_max, double arclen_inc
                , size_t successful_steps_max) {
    process_initialization(param_start,param_end,ds0);

    Matrix_t Y, Ynorm;
    Vector_t y = zeros_iteration<ParMethod>(x0,param_start);
    Vector_t ynorm;

    process_sub_iteration(y,Y,Ynorm,param_end,if_printiter);
    while (exitflag == EXITFLAG::OK) {
        ++processiter;
        y = method.template process<ParMethod>(predictor.predictor,y,_ds);
        
        if (!is_correct_solution(y,Y,arclen_min,arclen_inc)) continue;

        process_sub_iteration(y,Y,Ynorm,param_end,if_printiter);
        step_update(arclen_min,arclen_max,arclen_inc,successful_steps_max);
    }
    process_end_message();

}

template <typename Method, typename Pred>
void Continuation<Method,Pred>::process_sub_iteration(const Vector_t& y, Matrix_t& Y, Matrix_t& Ynorm, double param_end, bool if_printiter) {
    Vector_t ynorm = calculate_response_norm(y,method.ndof);
    
    Y.push_back(y);
    Ynorm.push_back(ynorm);

    predictor.calc_predictor(y,_ds);
    process_exitflag(y.last(), param_end);
    printiter(ynorm,y.last(),if_printiter);
}

template <typename Method, typename Pred>
void Continuation<Method,Pred>::process_initialization(double param_start, double param_end, double ds0) {
    _ds = ds0;
    processiter = 0;
    method.global_iter = 0;
    predictor.process_initialization();
    _successful_steps = 0;
    _count_not_correct_solution = 0;
}



template <typename Method, typename Pred>
template <typename ParMethod>
Vector_t Continuation<Method,Pred>::zeros_iteration(const Vector_t& x0, double param_start){
    Vector_t y0(x0.size()+1);
    Slice sl(y0.begin(),y0.end()-1,1);
    sl = x0;
    y0[x0.size()] = param_start;
    return method.template process<
        std::conditional_t<ParMethod::do_correction, ParMethod, typename corrector_methods::BasicCorrector>
    >(y0);
}

template <typename Method, typename Pred>
Vector_t Continuation<Method,Pred>::calculate_response_norm(const Vector_t& y, size_t ndof) {
    Vector_t ynorm(ndof);
    auto it = ynorm.begin();
    Slice<decltype(y.begin()),decltype(y.end())> yk;
    
    for (size_t k = 0; k < ndof; ++k, ++it) {
        yk.new_slice(y.begin()+k,y.end()-1,ndof);
        *it = norm(yk);
    }
    return ynorm;
}

template <typename Method, typename Pred>
void Continuation<Method,Pred>::process_exitflag(double param, double param_end) {
    if ((param > param_end) || param < 0) {
        exitflag = EXITFLAG::PARAM_END;
        return;
    }
    exitflag = EXITFLAG::OK;
}

template <typename Method, typename Pred>
void Continuation<Method,Pred>::printiter(const Vector_t& ynorm, double param, bool if_printiter) {
    if (if_printiter)
    std::cout 
            //   << "Iter: " << processiter
            //   << ", corr iter: " << method.iter
            //   << ", corr extflg: " << method.exitflag
            //   << ", {|y|}: [ " << ynorm << " ]"
            //   << ", param: " << param
              << ynorm << ' ' << param
              << std::endl;

}

template <typename Method, typename Pred>
void Continuation<Method,Pred>::process_end_message() {
    std::cout << "# Finished with exitflag: " << exitflag
              << ". Total number of corrector iterations: " << method.global_iter
              << "." << std::endl;
}

template <typename Method, typename Pred>
void Continuation<Method,Pred>::step_update(double arclen_min, double arclen_max, double arclen_inc
                                            , size_t successful_steps_max) {
    increase_step(arclen_max,arclen_inc,successful_steps_max);
}

template <typename Method, typename Pred>
bool Continuation<Method,Pred>::decrease_step(double arclen_min, double arclen_inc) {
    _successful_steps = 0;

    if (_ds  <= arclen_min)
        return false;
     
    _ds /= arclen_inc;
    std::cout << "# Arc-length step is decreased! New step size: "
              << _ds << std::endl;
    return true;
}

template <typename Method, typename Pred>
bool Continuation<Method,Pred>::increase_step(double arclen_max, double arclen_inc, size_t successful_steps_max) {
    if (_ds >= arclen_max)
        return false;
    if (_successful_steps < successful_steps_max) {
        ++_successful_steps;
        return false;
    }

    _ds *= arclen_inc;
    _successful_steps = 0;
    std::cout << "# Arc-length step is increased! New step size: "
              << _ds << std::endl;
    return true;

}

template <typename Method, typename Pred>
bool Continuation<Method,Pred>::is_correct_solution(const Vector_t& y, const Matrix_t& Y
            , double arclen_min, double arclen_inc)
{
    if (processiter < 2) return true;

    int flag = EXITFLAG::OK;
    if (method.exitflag == EXITFLAG::MAX_ITER) {
        flag = EXITFLAG::MAX_ITER;
    } else if ((norm(y-Y[Y.size()-2]) < method._epsx)) {
        flag = EXITFLAG::PREV_POINT;
    }

    if (flag == EXITFLAG::OK) {
        _count_not_correct_solution = 0;
        return true;
    }
    

    std::cout << "# Solution not found! FLAG " << flag << std::endl;
    ++_count_not_correct_solution;
    double ds_old = _ds;
    if (decrease_step(arclen_min, arclen_inc)) 
        predictor.resize_predictor(_ds/ds_old);


    return false;
}


} // namespace math