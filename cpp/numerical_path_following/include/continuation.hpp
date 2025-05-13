#pragma once

// std headers

// local headers
#include "linalg.hpp"
#include "numerical_continuation_base.hpp"

#include "corrector.hpp"
#include "basic_system.hpp"
#include "parametrization.hpp"
// other headers


namespace npath {

template <typename Corrector, typename Predictor>
class Continuation {

public:
    
    Corrector corrector;
    Predictor predictor;
private:
    double _ds;
    int exitflag;
    size_t processiter;
    size_t processiter_total;
    size_t _successful_steps;
    size_t _count_not_correct_solution;
    size_t _count_point_step_reduction;
    size_t max_not_correct_solution = 3;

public:

    Continuation() = default;
    Continuation(const Corrector& corrector, const Predictor& predictor);

    void process(math::vector_t<double,2>& Y,const math::vector<double>& y0, double param_start, double param_end, double ds0
                              , bool if_printiter, double arclen_min=1e-4, double arclen_max=1.0, double arclen_inc=2.0
                              , size_t successful_steps_max=10, size_t min_corrector_iters=4);
    
    void process_sub_iteration(const math::vector<double>& y
                                   , math::vector_t<double,2>& Y
                                   , math::vector_t<double,2>& Ynorm
                                   , math::vector<double>& ynorm
                                                , double param_end
                                                , bool if_printiter);
    
    
    void process_initialization(double param_start
                              , double param_end
                              , double ds0
                              , math::vector<double>& ynorm);

    
    void process_exitflag(double param, double param_end);
    void printiter(const math::vector<double>& ynorm, double param, const math::vector<double>& y, bool if_printiter);
    void process_end_message();
    void step_update(double arclen_min, double arclen_max, double arclen_inc
                    ,size_t successful_steps_max, size_t min_corrector_iters);
    bool decrease_step(double arclen_min, double arclen_inc);
    bool increase_step(double arclen_max, double arclen_inc, size_t successful_steps_max, size_t min_corrector_iters);
    bool is_correct_solution(math::vector<double>& y, math::vector_t<double,2>& Y, math::vector_t<double,2>& Ynorm, double arclen_min, double arclen_inc);
};


template <typename Corrector, typename Predictor>
Continuation<Corrector,Predictor>::Continuation(const Corrector& corrector, const Predictor& predictor)
        : corrector(corrector), predictor(predictor) {}



template <typename Corrector, typename Predictor>
void Continuation<Corrector,Predictor>::process(math::vector_t<double,2>& Y, const math::vector<double>& y0, double param_start, double param_end, double ds0
                , bool if_printiter, double arclen_min, double arclen_max, double arclen_inc
                , size_t successful_steps_max, size_t min_corrector_iters) {
    math::vector<double> ynorm;
    math::vector_t<double,2> Ynorm;
    math::vector<double> y = y0;
    process_initialization(param_start,param_end,ds0,ynorm);


    process_sub_iteration(y,Y,Ynorm,ynorm,param_end,if_printiter);
    while (exitflag == EXITFLAG::OK) {
        ++processiter_total;
        
        corrector.process(y,predictor.predictor,Y.last(),_ds);

        if (!is_correct_solution(y,Y,Ynorm,arclen_min,arclen_inc)) continue;

        ++processiter;
        
        step_update(arclen_min,arclen_max,arclen_inc,successful_steps_max, min_corrector_iters);
        process_sub_iteration(y,Y,Ynorm,ynorm,param_end,if_printiter);
    }
    process_end_message();
}

template <typename Corrector, typename Predictor>
void Continuation<Corrector,Predictor>::process_sub_iteration(const math::vector<double>& y, math::vector_t<double,2>& Y, math::vector_t<double,2>& Ynorm, math::vector<double>& ynorm, double param_end, bool if_printiter) {
    
    corrector.calculate_response_norm(y,ynorm);
    
    Y.push_back(y);
    Ynorm.push_back(ynorm);
    
    predictor.calc_predictor(_ds, Y);

    process_exitflag(y.last(), param_end);
    printiter(ynorm,y.last(), y, if_printiter);
}


template <typename Corrector, typename Predictor>
void Continuation<Corrector,Predictor>::process_initialization(double param_start
                                                             , double param_end
                                                             , double ds0
                                                             , math::vector<double>& ynorm) {
    _ds = ds0;
    processiter = 0;
    processiter_total = 0;
    corrector.continuation_initialization(ynorm);
    predictor.process_initialization();
    _successful_steps = 0;
    _count_not_correct_solution = 0;
    _count_point_step_reduction = 0;
}


template <typename Corrector, typename Predictor>
void Continuation<Corrector,Predictor>::process_exitflag(double param, double param_end) {
    if ((param > param_end) || param < 0) {
        exitflag = EXITFLAG::PARAM_END;
        return;
    }
    exitflag = EXITFLAG::OK;
}

template <typename Corrector, typename Predictor>
void Continuation<Corrector,Predictor>::printiter(const math::vector<double>& ynorm, double param, const math::vector<double>& y, bool if_printiter) {
    if (if_printiter)
    std::cout 
              << "\n#Iter: " << processiter
              << ", corr iter: " << corrector.iter
              << ", corr extflg: " << corrector.exitflag
              << ", {|y|}: [ " << ynorm << " ]"
              << ", param: " << param
              << ", step: " << _ds
            //   << "\n# predictor " << predictor.predictor
            //   << "\n# y " << y
              << '\n'
              << ynorm << ' ' << param
              << std::endl;

}

template <typename Corrector, typename Predictor>
void Continuation<Corrector,Predictor>::process_end_message() {
    std::cout << "# Finished with exitflag: " << exitflag
              << ". Number of continuation iterations: " << processiter
              << ". Total number of corrector iterations: " << corrector.global_iter
              << ". Total number of continuation iterations: " << processiter_total
              << ".\n";
}

template <typename Corrector, typename Predictor>
void Continuation<Corrector,Predictor>::step_update(double arclen_min, double arclen_max, double arclen_inc
                                            , size_t successful_steps_max, size_t min_corrector_iters) {
    increase_step(arclen_max,arclen_inc,successful_steps_max, min_corrector_iters);
}

template <typename Corrector, typename Predictor>
bool Continuation<Corrector,Predictor>::decrease_step(double arclen_min, double arclen_inc) {
    _successful_steps = 0;

    if (_ds  <= arclen_min)
        return false;
     
    _ds /= arclen_inc;
    std::cout << "# Arc-length step is decreased! New step size: "
              << _ds << std::endl;
    return true;
}

template <typename Corrector, typename Predictor>
bool Continuation<Corrector,Predictor>::increase_step(double arclen_max, double arclen_inc, size_t successful_steps_max, size_t min_corrector_iters) {
    if (_ds >= arclen_max || corrector.iter > min_corrector_iters)
        return false;
    if (_successful_steps < successful_steps_max) {
        ++_successful_steps;
        return false;
    }

    _ds *= arclen_inc;
    _successful_steps = 0;
    std::cout << "# Arc-length step is increased! New step size: "
              << _ds << '\n';
    return true;

}

template <typename Corrector, typename Predictor>
bool Continuation<Corrector,Predictor>::is_correct_solution(math::vector<double>& y, math::vector_t<double,2>& Y, math::vector_t<double,2>& Ynorm
            , double arclen_min, double arclen_inc)
{
    if (processiter < 2) return true;

    int flag = EXITFLAG::OK;
    if (corrector.exitflag == EXITFLAG::MAX_ITER) {
        flag = EXITFLAG::MAX_ITER;
    } else if ((norm(y-Y[Y.size()-2]) < corrector._epsx)) {
        flag = EXITFLAG::PREV_POINT;
    }

    if (flag == EXITFLAG::OK) {
        _count_not_correct_solution = 0;
        _count_point_step_reduction = 0;
        return true;
    }
    std::cout << "# Solution not found! FLAG " << flag << '\n';
    ++_count_not_correct_solution;
    
    if (_count_not_correct_solution == max_not_correct_solution) {
        _count_not_correct_solution = 0;
        if (_count_point_step_reduction > 1)
            _ds *= pow(arclen_inc,_count_point_step_reduction-1);
        std::cout << "# Step back! New step length " << _ds << std::endl;
        _count_point_step_reduction = 0;
        Y.erase(Y.end()-1);
        Ynorm.erase(Ynorm.end()-1);
        corrector.step_back(_ds);
        predictor.step_back();
        predictor.calc_predictor(_ds, Y);
        y = *(Y.end()-1);

    } else {
        double ds_old = _ds;
        if (decrease_step(arclen_min, arclen_inc)) {
            predictor.resize_predictor(_ds,ds_old);
            ++_count_point_step_reduction;
        }
        corrector.solution_not_found();
    }

        

    return false;
}

template <typename Corrector, typename Predictor>
Continuation(const Corrector&, const Predictor&) -> Continuation<Corrector,Predictor>;


} // namespace math