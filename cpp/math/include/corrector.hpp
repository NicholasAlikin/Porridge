/*
Corrector algorithms:
    o default: Newton-Raphson algorithm
    o normal-flow: Gauss-Newton algorithm
                   using Moore-Penrose inverse
                   (pseudoinverse) matrix
    o Psevdo-arc-lenght algorithms:
        - sphere scheme,
        - orthogonal scheme.
*/
#pragma once

#include "linalg.hpp"
#include "dft.hpp"

namespace math { 

using Vector_t = vector_t<double,1>;
using Matrix_t = vector_t<double,2>;


class Corrector {

public:
    static constexpr int EXTFLG_OK = 0;
    static constexpr int EXTFLG_VAR_FUN = 1;
    static constexpr int EXTFLG_MAX_ITER = 2;
    
    

    size_t iter = 0, global_iter = 0;
    int exitflag = 0;
    double _dy_norm = 0, _fun_norm = 0;
    double _epsx = 1e-5, _epsf = 1e-5;
    size_t _max_iter = 20;

    Matrix_t mass;
    Matrix_t damp;
    Matrix_t stif;

    size_t ndof;
    DFT dft;
    

    Corrector() = default;
    Corrector(const Matrix_t& mass
            , const Matrix_t& damp
            , const Matrix_t& stif
            , const DFT& dft);
    Corrector(const Matrix_t& mass
            , const Matrix_t& damp
            , const Matrix_t& stif
            , size_t N, size_t H);

    ~Corrector() {}

private:

    void process_initialization();
    
    void process_callback();
    void process_exitflag(const Vector_t& dy, const Vector_t& fun);
    void printiter();

    void correction_sub_iteration(const Vector_t& y, const Vector_t& dy, const Vector_t& fun, const Matrix_t& jac);
    void correction_system_response(const Vector_t& fun, const Matrix_t& jac, const Vector_t& y);
    

public:

    // based on corrector algorithm
    Vector_t correction(const Vector_t& predictor); // virtual ???
    Vector_t continuation(const Vector_t& predictor, const Vector_t& previous); // virtual ???

    Vector_t process_iteration(const Matrix_t& jac, const Vector_t& fun);
    
    
    void correction_total_increment(const Vector_t& y, const Vector_t& dy);
    
};


} // namespace math