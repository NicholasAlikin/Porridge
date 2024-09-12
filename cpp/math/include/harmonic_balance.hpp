#pragma once

#include "corrector.hpp"
#include <functional>
#include "dft.hpp"

namespace math {


class HBM: public Corrector {
private:
public:
    Matrix_t Mass;
    Matrix_t Damp;
    Matrix_t Stif;
    std::function<Vector_t(Vector_t,double,DFT)> funcnl;
    std::function<Vector_t(double,DFT)> funcex;

    HBM() = default;
    
    HBM(const Matrix_t& mass
      , const Matrix_t& damp
      , const Matrix_t& stif
      , const DFT& dft
      , const std::function<Vector_t(const Vector_t&,double,const DFT&)>& funcnl
      , const std::function<Vector_t(double,const DFT&)>& funcex);

    void correction_system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y) override;
    
    
    Matrix_t linear_system_dynamic_reaction(double freq);

    Vector_t correction_system(const Vector_t& x, const Matrix_t& L, const Vector_t& fnl, const Vector_t& fex);
    Matrix_t corretion_system_jac(Vector_t& x, double freq, const Matrix_t& L, const Vector_t& fnl);
};

} // namespace math