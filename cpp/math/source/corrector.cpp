#include "corrector.hpp"

namespace math {

Corrector::Corrector(const Matrix_t& mass
                    ,const Matrix_t& damp
                    ,const Matrix_t& stif
                    ,const DFT& dft)
        : mass(mass), damp(damp), stif(stif)
        , dft(dft) {

}

Corrector::Corrector(const Matrix_t& mass
                    ,const Matrix_t& damp
                    ,const Matrix_t& stif
                    ,size_t N, size_t H)
        : mass(mass), damp(damp), stif(stif)
        , dft(DFT(N,H,mass.size())) {

}

Vector_t Corrector::correction(const Vector_t& predictor) {
    process_initialization();
    Vector_t y = predictor;
    Vector_t dy;
    Vector_t fun = zeros<double>(y.size());
    Matrix_t jac = zeros<double>(y.size(),y.size());

    correction_sub_iteration(y,dy,fun,jac);

    while (exitflag == Corrector::EXTFLG_OK) {
        dy = process_iteration(fun,jac, dy);
        correction_total_increment(y,dy);
        ++iter; ++global_iter;

        correction_sub_iteration(y,dy,fun,jac);
    }

    return y;
}

void Corrector::process_initialization() {
    this->iter = 0;
}

void Corrector::process_exitflag(const Vector_t& dy, const Vector_t& fun) {
    _dy_norm = norm(dy);
    _fun_norm = norm(fun);
    if (_dy_norm < _epsx && _fun_norm < _epsf) {
        exitflag = Corrector::EXTFLG_VAR_FUN;
    } else if (iter > _max_iter) {
        exitflag = Corrector::EXTFLG_MAX_ITER;
    } else {
        exitflag = Corrector::EXTFLG_OK;
    }
}

void Corrector::printiter() {
    std::cout << "Iter: " << iter
              << ", |f|: " << _fun_norm
              << ", |dy|: " << _dy_norm
              << '\n';
}

Vector_t math::Corrector::process_iteration(const Matrix_t &jac, const Vector_t &fun)
{
    return solve(jac,fun);
}


void Corrector::correction_sub_iteration(const Vector_t& y, const Vector_t& dy, const Vector_t& fun, const Matrix_t& jac) {
    correction_system_response(fun,jac, y);
    process_exitflag(dy,fun);
    printiter();
}




} // namespace math