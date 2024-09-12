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
        dy = process_iteration(jac,fun);
        correction_total_increment(y,dy);
        ++iter; ++global_iter;

        correction_sub_iteration(y,dy,fun,jac);
    }

    return y;
}

Vector_t Corrector::continuation(const Vector_t& predictor, const Vector_t& previous) {
    process_initialization();
    Vector_t y = predictor + previous;
    Vector_t dy; // zeros initializetion is not neaded beacuse norm of empty vector = 0
    Vector_t Dy = predictor;
    Vector_t fun = zeros<double>(y.size());
    Matrix_t jac = zeros<double>(y.size(),y.size());

    continuation_sub_iteration(y,dy,Dy,fun,jac);

    while (exitflag == Corrector::EXTFLG_OK) {
        dy = process_iteration(jac,fun);
        continuation_total_increment(Dy,dy,predictor);
        y = previous + Dy;
        ++iter; ++global_iter;

        continuation_sub_iteration(y,dy,Dy,fun,jac);
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
              << std::endl;
}

Vector_t Corrector::process_iteration(const Matrix_t &jac, const Vector_t& fun)
{
    return solve(jac,fun);
}


void Corrector::correction_sub_iteration(const Vector_t& y, const Vector_t& dy, Vector_t& fun, Matrix_t& jac) {
    correction_system_response(fun,jac, y);
    process_exitflag(dy,fun);
    printiter();
}
void Corrector::continuation_sub_iteration(const Vector_t& y, const Vector_t& dy, const Vector_t& Dy, Vector_t& fun, Matrix_t& jac) {
    continuation_system_response(fun,jac, y,Dy);
    process_exitflag(dy,fun);
    printiter();
}

void Corrector::correction_total_increment(Vector_t& y, const Vector_t& dy) {
    Slice sl(y.begin(),y.end()-1);
    sl += dy;
}

void Corrector::continuation_total_increment(Vector_t& Dy, const Vector_t& dy, const Vector_t& predictor) {
    Dy += dy;
}

void Corrector::continuation_system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t& Dy) {
    correction_system_response(fun,jac, y);
}

void NormalFlow::correction_total_increment(Vector_t& y, const Vector_t& dy) {
    y += dy;
}

Vector_t NormalFlow::process_iteration(const Matrix_t& jac, const Vector_t& fun) {
    return dot(pinv(jac),fun);
}

} // namespace math