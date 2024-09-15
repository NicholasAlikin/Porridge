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

namespace corrector_methods {
struct BasicCorrector;
struct NormalFlowCorrector;
struct ArcLengthCorrector;
} // namespace corrector_methods


class Corrector {
public:
    static constexpr int EXTFLG_OK = 0;
    static constexpr int EXTFLG_VAR_FUN = 1;
    static constexpr int EXTFLG_MAX_ITER = 2;



    mutable size_t iter = 0, global_iter = 0;
    mutable int exitflag = 0;
    mutable double _dy_norm = 0, _fun_norm = 0;
    double _epsx = 1e-5, _epsf = 1e-5;
    double _dx = 1e-5;
    size_t _max_iter = 10;

    Matrix_t mass;
    Matrix_t damp;
    Matrix_t stif;

    size_t ndof;
    DFT dft;
    


    Corrector(const Matrix_t& mass
            , const Matrix_t& damp
            , const Matrix_t& stif
            , const DFT& dft);
    Corrector(const Matrix_t& mass
            , const Matrix_t& damp
            , const Matrix_t& stif
            , size_t N, size_t H);

    ~Corrector() = default;

public:

    void process_initialization() const;
    
    void process_callback()  const;
    void process_exitflag(const Vector_t& dy, const Vector_t& fun) const;
    void printiter() const;

    template <typename Method>
    void process_sub_iteration(const Vector_t& y, const Vector_t& dy, Vector_t& fun, Matrix_t& jac) const;
    template <typename Method>
    void process_sub_iteration(const Vector_t& y, const Vector_t& dy, const Vector_t& Dy, Vector_t& fun, Matrix_t& jac, double ds) const;
    
    virtual void system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y) const = 0;
    virtual void system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t& Dy) const = 0;
    virtual void system_response_ext(Vector_t& fun, Matrix_t& jac, const Vector_t& y) const = 0;
    virtual void system_response_ext(Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t& Dy) const = 0;
    

public:

    // based on corrector algorithm
    template <typename Method = corrector_methods::BasicCorrector>
    Vector_t process(const Vector_t& predictor) const;
    template <typename Method = corrector_methods::BasicCorrector>
    Vector_t process(const Vector_t& predictor, const Vector_t& previous, double ds) const;
    
};

template <typename Method>
Vector_t Corrector::process(const Vector_t& inital) const {
    process_initialization();
    Vector_t y = inital;
    Vector_t dy;
    Vector_t fun = Method::get_init_fun(y.size());
    Matrix_t jac = Method::get_init_jac(y.size());

    process_sub_iteration<Method>(y,dy,fun,jac);

    while (exitflag == Corrector::EXTFLG_OK) {
        dy = Method::process_iteration(jac,fun);
        Method::process_total_increment(y,dy);
        ++iter; ++global_iter;

        process_sub_iteration<Method>(y,dy,fun,jac);
    }

    return y;
}

template <typename Method>
Vector_t Corrector::process(const Vector_t& predictor, const Vector_t& previous, double ds) const {
    process_initialization();
    Vector_t y = predictor + previous;
    Vector_t dy; // zeros initializetion is not neaded beacuse norm of empty vector = 0
    Vector_t Dy = predictor;
    Vector_t fun = Method::get_init_fun(y.size());
    Matrix_t jac = Method::get_init_jac(y.size());
    
    process_sub_iteration<Method>(y,dy,Dy,fun,jac,ds);
    while (exitflag == Corrector::EXTFLG_OK) {
        dy = Method::process_iteration(jac,fun);
        Method::process_total_increment(Dy,dy,predictor,ds);
        y = previous + Dy;
        ++iter; ++global_iter;
        process_sub_iteration<Method>(y,dy,Dy,fun,jac,ds);
    }

    return y;
}

template <typename Method>
void Corrector::process_sub_iteration(const Vector_t& y, const Vector_t& dy, Vector_t& fun, Matrix_t& jac) const {
    Method::system_response(*this, fun,jac, y);
    process_exitflag(dy,fun);
    printiter();
}
template <typename Method>
void Corrector::process_sub_iteration(const Vector_t& y, const Vector_t& dy, const Vector_t& Dy, Vector_t& fun, Matrix_t& jac, double ds) const {
    Method::system_response(*this, fun,jac, y,Dy,ds);
    process_exitflag(dy,fun);
    printiter();
}

namespace corrector_methods {

struct BasicCorrector {
    static const bool do_correction = true;
    static const bool do_continuation = false;

    static Vector_t get_init_fun(size_t y_size);
    static Matrix_t get_init_jac(size_t y_size);

    static Vector_t process_iteration(const Matrix_t& jac, const Vector_t& fun);

    static void process_total_increment(Vector_t& y, const Vector_t& dy);
    static void process_total_increment(Vector_t& Dy, const Vector_t& dy, const Vector_t& predictor, double ds);
    
    static void system_response(const Corrector& corrector, Vector_t& fun, Matrix_t& jac, const Vector_t& y);
    static void system_response(const Corrector& corrector, Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t& Dy, double ds) = delete;

};


struct NormalFlowCorrector: BasicCorrector {
    static const bool do_continuation = true;
    
    static Matrix_t get_init_jac(size_t y_size);
    static Vector_t process_iteration(const Matrix_t& jac, const Vector_t& fun);

    static void process_total_increment(Vector_t& y, const Vector_t& dy);
    static void process_total_increment(Vector_t& Dy, const Vector_t& dy, const Vector_t& predictor, double ds);

    static void system_response(const Corrector& corrector, Vector_t& fun, Matrix_t& jac, const Vector_t& y);
    static void system_response(const Corrector& corrector, Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t& Dy, double ds);
};

struct ArcLengthCorrector: BasicCorrector {
    static const bool do_correction = false;
    static const bool do_continuation = true;

    static Vector_t get_init_fun(size_t y_size);
    static Matrix_t get_init_jac(size_t y_size);

    static void process_total_increment(Vector_t& y, const Vector_t& dy) = delete;
    static void process_total_increment(Vector_t& Dy, const Vector_t& dy, const Vector_t& predictor, double ds);

    static void system_response(const Corrector& corrector, Vector_t& fun, Matrix_t& jac, const Vector_t& y) = delete;
    static void system_response(const Corrector& corrector, Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t& Dy, double ds);


};

} // namespace corrector_methods 

// class ArcLength: public Corrector {
// public:

// };

} // namespace math