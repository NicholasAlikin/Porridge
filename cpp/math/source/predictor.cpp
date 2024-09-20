#include "predictor.hpp"

namespace math {

Predictor::Predictor(size_t size): predictor(zeros<double>(size)) {}

void Predictor::calc_predictor(double ds) {
    predictor[predictor.size()-1] = ds;
}

void Predictor::process_initialization() {
    auto it = predictor.begin(), end = predictor.end();
    while (it != end) {
        *it = 0;
        ++it;
    }
}

Secant::Secant(size_t size): Predictor(size), previous(zeros<double>(size)) {}

void Secant::calc_predictor(const Vector_t& y, double ds) {
    Vector_t prev = y - std::move(previous);
    
    double prev_norm = norm(prev);
    predictor = std::move(prev);
    predictor *= ds/prev_norm;
    
    previous = y;

    // predictor.last() = ds;
}

void Secant::process_initialization() {
    this->Predictor::process_initialization();
    auto it = previous.begin(), end = previous.end();
    while (it != end) {
        *it = 0;
        ++it;
    }
}

} // namespace math 