#include "corrector.hpp"

namespace math {

Corrector::Corrector(const Matrix_t& mass
                    ,const Matrix_t& damp
                    ,const Matrix_t& stif
                    ,const DFT& dft)
        : mass(mass), damp(damp), stif(stif)
        , ndof(dft.ndof), dft(dft) {

}

Corrector::Corrector(const Matrix_t& mass
                    ,const Matrix_t& damp
                    ,const Matrix_t& stif
                    ,size_t N, size_t H)
        : mass(mass), damp(damp), stif(stif)
        , ndof(mass.size()), dft(DFT(N,H,mass.size())) {

}





void Corrector::process_initialization() const {
    iter = 0;
}

void Corrector::process_exitflag(const Vector_t& dy, const Vector_t& fun) const {
    _dy_norm = norm(dy);
    _fun_norm = norm(fun);
    if (_dy_norm < _epsx && _fun_norm < _epsf) {
        exitflag = Corrector::EXTFLG_VAR_FUN;
    } else if (iter >= _max_iter) {
        exitflag = Corrector::EXTFLG_MAX_ITER;
    } else {
        exitflag = Corrector::EXTFLG_OK;
    }
}

void Corrector::printiter() const {
    // std::cout << "Iter: " << iter
    //           << ", |f|: " << _fun_norm
    //           << ", |dy|: " << _dy_norm
    //           << std::endl;
}






namespace corrector_methods {
/*------------
BasicCorrector
------------*/

Vector_t BasicCorrector::get_init_fun(size_t y_size) {
    Vector_t fun = zeros<double>(y_size-1);
    return fun;
}

Matrix_t BasicCorrector::get_init_jac(size_t y_size) {
    Matrix_t jac = zeros<double>(y_size-1,y_size-1);
    return jac;
}

Vector_t BasicCorrector::process_iteration(const Matrix_t& jac, const Vector_t& fun) {
    return solve(jac,fun);
}

void BasicCorrector::process_total_increment(Vector_t& y, const Vector_t& dy) {
    Slice sl(y.begin(),y.end()-1);
    sl += dy;
}
void BasicCorrector::process_total_increment(Vector_t& Dy, const Vector_t& dy, const Vector_t& predictor, double ds) {
    BasicCorrector::process_total_increment(Dy,dy);
}

void BasicCorrector::system_response(const Corrector& corrector, Vector_t& fun, Matrix_t& jac, const Vector_t& y) {
    corrector.system_response(fun,jac, y);
}


/*------------
NormalFlowCorrector
------------*/
Matrix_t NormalFlowCorrector::get_init_jac(size_t y_size) {
    Matrix_t jac = zeros<double>(y_size-1,y_size);
    return jac;
}

Vector_t NormalFlowCorrector::process_iteration(const Matrix_t& jac, const Vector_t& fun) {
    return psolve(jac,fun);
}

void NormalFlowCorrector::process_total_increment(Vector_t& y, const Vector_t& dy) {
    y += dy;
}
void NormalFlowCorrector::process_total_increment(Vector_t& Dy, const Vector_t& dy, const Vector_t& predictor, double ds) {
    Dy += dy;
    Dy *= ds/norm(Dy);
}

void NormalFlowCorrector::system_response(const Corrector& corrector, Vector_t& fun, Matrix_t& jac, const Vector_t& y) {
    corrector.system_response_ext(fun,jac, y);
}

void NormalFlowCorrector::system_response(const Corrector& corrector, Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t&, double) {
    NormalFlowCorrector::system_response(corrector, fun,jac, y);
}

/*------------
ArcLengthCorrector
------------*/

Vector_t ArcLengthCorrector::get_init_fun(size_t y_size) {
    Vector_t fun = zeros<double>(y_size);
    return fun;
}

Matrix_t ArcLengthCorrector::get_init_jac(size_t y_size) {
    Matrix_t jac = zeros<double>(y_size,y_size);
    return jac;
}


void ArcLengthCorrector::process_total_increment(Vector_t& Dy, const Vector_t& dy, const Vector_t& predictor, double ds) {
    Dy += dy;
}

void ArcLengthCorrector::system_response(const Corrector& corrector, Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t& Dy, double ds) {
    corrector.system_response_ext(fun,jac, y);
    fun.last() = -(dot(Dy,Dy) - ds*ds);
    jac.last() = 2*Dy;
}



} // namespace corrector_methods 
} // namespace math