#include "predictor.hpp"

namespace npath {

Predictor::Predictor(size_t size): _sz(size) {}
// Predictor::Predictor(const Predictor& other): _sz(other._sz) {}

void Predictor::calc_predictor(double ds, const math::vector_t<double,2>&) {
    predictor.last() = ds;
}

void Predictor::process_initialization() {
    predictor = math::zeros<double>(_sz);
}

void Predictor::resize_predictor(double new_step, double) {
    predictor.last() = new_step;
}


void Secant::calc_predictor(double ds, const math::vector_t<double,2>& Y) {
    
    auto itpre = predictor.begin(), endpre = predictor.end();
    auto itcur = Y.last().begin();

    // if previous point was calculated
    double prenorm = 0;
    
    if (Y.size() > 1) {
        auto itprv = (*(Y.end()-2)).begin(); // previuos
        while (itpre < endpre) {
            *itpre = (*itcur) - (*itprv);         // only once value will be loaded
            prenorm += (*itpre) * (*itpre); // from the `itpre` address
            ++itpre;
            ++itcur;
            ++itprv;
        }
    } else {
        // previous = 0
        while (itpre < endpre) {
            *itpre = (*itcur);         // only once value will be loaded
            prenorm += (*itpre) * (*itpre); // from the `itpre` address
            ++itpre;
            ++itcur;
        }
    }

    predictor *= ds/sqrt(prenorm);
}

void Secant::resize_predictor(double new_step, double old_step) {
    predictor *= new_step/old_step;
}


void Polynomical::process_initialization() {
    this->Predictor::process_initialization();
    steps_history = math::vector<double>();
    coefs = math::vector_t<double,2>(3);
}


void Polynomical::calc_predictor(double ds, const math::vector_t<double,2>& Y) {
    // no saved steps
    
    
    if (steps_history.size() == 0) {
        steps_history.push_back(norm(*(Y.end()-2) - *(Y.end()-3)));
    }
    steps_history.push_back(norm(*(Y.end()-1)  - *(Y.end()-2) ));
       
    auto last_step = steps_history.end();
    double ds_prv = *(last_step-1), ds_prv2 = *(last_step-2);
    coefs[0] = *(Y.end()-3);
    double tmp = 1.0/ (ds_prv*ds_prv2*(ds_prv+ds_prv2));
    
    coefs[1] =  ( tmp*(ds_prv+ds_prv2)*(ds_prv+ds_prv2) )*(*(Y.end()-2)) 
               -( tmp* ds_prv*(ds_prv+2*ds_prv2)        )*(*(Y.end()-3))
               -( tmp* ds_prv2*ds_prv2                  )*(*(Y.end()-1));
    
    coefs[2] = -( tmp*(ds_prv+ds_prv2)                  )*(*(Y.end()-2))   
               +( tmp* ds_prv                           )*(*(Y.end()-3))
               +( tmp* ds_prv2                          )*(*(Y.end()-1));
    
    predictor = coefs[0] + coefs[1]*(ds_prv+ds_prv2+ds) + coefs[1]*(ds_prv+ds_prv2+ds)*(ds_prv+ds_prv2+ds); 
}

void Polynomical::resize_predictor(double new_step, double) {
    auto last_step = steps_history.end();
    double ds_prv = *(last_step-1), ds_prv2 = *(last_step-2);

    predictor = coefs[0] + coefs[1]*(ds_prv+ds_prv2+new_step) + coefs[1]*(ds_prv+ds_prv2+new_step)*(ds_prv+ds_prv2+new_step);
}

void Polynomical::step_back() {
    steps_history.erase(steps_history.end()-1);
}


} // namespace npath 