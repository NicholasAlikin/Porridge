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
    double _dx = 1e-5;
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

    void correction_sub_iteration(const Vector_t& y, const Vector_t& dy, Vector_t& fun, Matrix_t& jac);
    void continuation_sub_iteration(const Vector_t& y, const Vector_t& dy, const Vector_t& Dy, Vector_t& fun, Matrix_t& jac);
    
    virtual void correction_system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y) {};
    virtual void continuation_system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t& Dy);
    

public:

    // based on corrector algorithm
    Vector_t correction(const Vector_t& predictor); // virtual ???
    Vector_t continuation(const Vector_t& predictor, const Vector_t& previous); // virtual ???

    virtual Vector_t process_iteration(const Matrix_t& jac, const Vector_t& fun);
    
    
    virtual void correction_total_increment(Vector_t& y, const Vector_t& dy);
    void continuation_total_increment(Vector_t& Dy, const Vector_t& dy, const Vector_t& predictor);
    
};


class NormalFlow: public Corrector {
public:
    void correction_total_increment(Vector_t& y, const Vector_t& dy) override;
    Vector_t process_iteration(const Matrix_t& jac, const Vector_t& fun) override;

};

// class ArcLength: public Corrector {
// public:

// };

} // namespace math