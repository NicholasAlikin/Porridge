#pragma once

#include "linalg.hpp"

namespace npath {

template <typename Pred>
struct PredictorTraits;



class Predictor {
public:
    math::vector<double> predictor; // predictor value
    size_t _sz; // predictor size

    Predictor() = default;
    // Predictor(const Predictor& other);
    Predictor(size_t size);

    void calc_predictor(double step, const math::vector_t<double,2>& Y);
    void process_initialization();
    void resize_predictor(double new_step, double old_step);
    void step_back() {}
    
};



class Secant: public Predictor {
public:
    using Predictor::Predictor;
    void calc_predictor(double ds, const math::vector_t<double,2>& Y);

    void resize_predictor(double new_step, double old_step);
};


class Polynomical: public Secant {
public:
    using Secant::Secant;
    
    
    math::vector<double> steps_history;
    math::vector_t<double,2> coefs;

    void process_initialization();
    void calc_predictor(double step, const math::vector_t<double,2>& Y);
    void resize_predictor(double new_step, double old_step);
    void step_back();
};


} // namespace npath