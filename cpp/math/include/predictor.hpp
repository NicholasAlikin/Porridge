#pragma once

#include "linalg.hpp"

namespace math {

using Vector_t = vector<double>;

class Predictor {
public:
    Vector_t predictor;
    Predictor() = default;
    Predictor(size_t size);

    void calc_predictor(double ds);
    void process_initialization();
    void resize_predictor(double ds);
    
};

class Secant: public Predictor {
public:
    Vector_t previous;
    Secant() = default;
    Secant(size_t size);

    void calc_predictor(const Vector_t& y, double ds);
    void process_initialization();
};

} // namespace math