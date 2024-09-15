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

    
    HBM(const Matrix_t& mass
      , const Matrix_t& damp
      , const Matrix_t& stif
      , const DFT& dft
      , const std::function<Vector_t(const Vector_t&,double,const DFT&)>& funcnl
      , const std::function<Vector_t(double,const DFT&)>& funcex);


    void system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y) const override;
    void system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t& Dy) const override;
    void system_response_ext(Vector_t& fun, Matrix_t& jac, const Vector_t& y) const override;
    void system_response_ext(Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t& Dy) const override;

    // void correction_system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y) override;
    // void continuation_system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t& Dy) override;

    
    Matrix_t linear_system_dynamic_reaction(double freq) const;

    Vector_t system(const Vector_t& x, const Matrix_t& L, const Vector_t& fnl, const Vector_t& fex) const;
    void system_jac(Matrix_t& jac, Vector_t& x, double freq, const Matrix_t& L, const Vector_t& fnl) const;
    void system_jac_freq(Vector_t& jac, const Vector_t& x, double freq, const Vector_t& fnl, const Vector_t& fex) const;
};

} // namespace math