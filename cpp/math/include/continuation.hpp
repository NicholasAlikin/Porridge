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



    Continuation() = default;
    Continuation(const Method& method, const Pred& predictor);

    template <typename ParMethod>
    requires ParMethod::do_continuation
    void process(const Vector_t& x0, double param_start, double param_end, double ds0);
    
    void process_initialization(double param_start, double param_end, double ds0);

    template <typename ParMethod>
    Vector_t zeros_iteration(const Vector_t& x0,double param_start);
    
    void process_exitflag(double param, double param_end);
    void printiter(const Vector_t& ynorm, double param);
    Vector_t calculate_response_norm(const Vector_t& y, size_t ndof);
    void process_end_message();
};


template <typename Method, typename Pred>
Continuation<Method,Pred>::Continuation(const Method& method, const Pred& predictor)
        : method(method), predictor(predictor)
{

}



template <typename Method, typename Pred>
template <typename ParMethod>
requires ParMethod::do_continuation
void Continuation<Method,Pred>::process(const Vector_t& x0, double param_start, double param_end, double ds0) {
    process_initialization(param_start,param_end,ds0);

    Matrix_t Y, Ynorm;
    Vector_t y = zeros_iteration<ParMethod>(x0,param_start);
    Vector_t ynorm = calculate_response_norm(y,method.ndof);
    Y.push_back(y);
    Ynorm.push_back(ynorm);

    predictor.calc_predictor(y,_ds);
    process_exitflag(y.last(), param_end);
    printiter(ynorm,y.last());
    while (exitflag == EXITFLAG::OK) {
        ++processiter;
        y = method.template process<ParMethod>(predictor.predictor,y,_ds);
        
        ynorm = calculate_response_norm(y,method.ndof);
        
        Y.push_back(y);
        Ynorm.push_back(ynorm);
        predictor.calc_predictor(y,_ds);
        
        process_exitflag(y.last(), param_end);
        printiter(ynorm,y.last());
    }
    process_end_message();

}

template <typename Method, typename Pred>
void Continuation<Method,Pred>::process_initialization(double param_start, double param_end, double ds0) {
    _ds = ds0;
    processiter = 0;
    method.global_iter = 0;
    predictor.process_initialization();
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
void Continuation<Method,Pred>::printiter(const Vector_t& ynorm, double param) {
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
    std::cout << "Finished with exitflag: " << exitflag
              << ". Total number of corrector iterations: " << method.global_iter
              << "." << std::endl;
}

} // namespace math