#pragma once
#include "linalg.hpp"

namespace math {

using Matrix_t = vector_t<double, 2>;
using Vector_t = vector_t<double, 1>;

class DFT {

public:
    size_t H = 0;
    size_t N = 0;
    size_t ndof = 1;

    // Matrix_t forward_basic;     // for ndof = 1
    Matrix_t forward_basic;
    Matrix_t backward_basic;
    Matrix_t derivative_basic;  // for ndof = 1
    Matrix_t derivative2_basic; // for ndof = 1

    // Matrix_t forward;
    // Matrix_t backward;
    // Matrix_t derivative;
    // Matrix_t derivative2;


    static Matrix_t calculate_forward_basic(size_t H, size_t N);
    static Matrix_t calculate_backward_basic(size_t H, size_t N);
    static Matrix_t calculate_derivative_basic(size_t H, size_t N);
    static Matrix_t calculate_derivative2_basic(size_t H, size_t N);
    static Matrix_t calculate_derivative2_basic(const Matrix_t& T_derivative_basic);


    static Matrix_t calculate_forward(size_t H, size_t N, size_t ndof);
    static Matrix_t calculate_backward(size_t H, size_t N, size_t ndof);
    static Matrix_t calculate_derivative(size_t H, size_t N, size_t ndof);


    DFT() = default;
    DFT(size_t H, size_t N, size_t ndof = 1);

private:
    Vector_t transform(const Matrix_t& transform_matrix, size_t res_size, const Vector_t& vec) const;
    void __derivative(Vector_t& res) const;
    void __derivative2(Vector_t& res) const;
public:
    Vector_t forward(const Vector_t& vec) const;
    Vector_t backward(const Vector_t& vec) const;
    Vector_t derivative(const Vector_t& vec) const;
    Vector_t derivative(Vector_t&& vec) const;
    Vector_t derivative2(const Vector_t& vec) const;
    Vector_t derivative2(Vector_t&& vec) const;
    

};

} // namespace math