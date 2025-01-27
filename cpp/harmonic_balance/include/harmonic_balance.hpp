#pragma once

#include <functional>

#include "basic_system.hpp"
#include "dft.hpp"

namespace npath {


class HBM: public BasicSystem {
private:
public:
    
    DFT dft;
    size_t ndof;
    math::vector_t<double,2> Mass;
    math::vector_t<double,2> Damp;
    math::vector_t<double,2> Stif;
    std::function<math::vector<double>(math::vector<double>,double,DFT)> funcnl;
    std::function<math::vector<double>(double,DFT)> funcex;

    double _dx = 1e-5;

    HBM(const math::vector_t<double,2>& mass
      , const math::vector_t<double,2>& damp
      , const math::vector_t<double,2>& stif
      , const DFT& dft
      , const std::function<math::vector<double>(const math::vector<double>&,double,const DFT&)>& funcnl);
    
    HBM(const math::vector_t<double,2>& mass
      , const math::vector_t<double,2>& damp
      , const math::vector_t<double,2>& stif
      , const DFT& dft
      , const std::function<math::vector<double>(const math::vector<double>&,double,const DFT&)>& funcnl
      , const std::function<math::vector<double>(double,const DFT&)>& funcex);

    void continuation_initialization(math::vector_t<double,2>& jac
                                     ,math::vector<double>& ynorm);

    void system_response(math::vector<double>& fun, math::vector_t<double,2>& jac, math::vector<double>& y);
    void system_response(math::vector<double>& fun, math::vector_t<double,2>& jac, math::vector<double>& y, const math::vector<double>& Dy);
    void system_response_extended(math::vector<double>& fun, math::vector_t<double,2>& jac, math::vector<double>& y);
    void system_response_extended(math::vector<double>& fun, math::vector_t<double,2>& jac, math::vector<double>& y, const math::vector<double>& Dy);


    void calculate_response_norm(const math::vector<double>& y
                                     , math::vector<double>& ynorm);

    size_t response_norm_size();  


    math::vector_t<double,2> linear_system_dynamic_reaction(double freq) const;

    math::vector<double> system(const math::vector<double>& x, const math::vector_t<double,2>& L, const math::vector<double>& fnl, const math::vector<double>& fex) const;
    void system_jac(math::vector_t<double,2>& jac, math::vector<double>& x, double freq, const math::vector_t<double,2>& L, const math::vector<double>& fnl) const;
    void system_jac_freq(math::vector<double>& jac, const math::vector<double>& x, double freq, const math::vector<double>& fnl, const math::vector<double>& fex) const;
};


} // namespace npath