#pragma once
#include "linalg.hpp"

namespace math {

using Matrix_t = vector_t<double, 2>;

class DFT {

public:
    size_t H = 0;
    size_t N = 0;
    size_t ndof = 1;

    Matrix_t forward;
    Matrix_t backward;
    Matrix_t derivative;
    Matrix_t derivative2;

    Matrix_t forward_basic;     // for ndof = 1
    Matrix_t derivative_basic;  // for ndof = 1
    Matrix_t derivative2_basic; // for ndof = 1

    static Matrix_t calculate_forward_basic(size_t N, size_t H);
    static Matrix_t calculate_backward_basic(size_t N, size_t H);
    static Matrix_t calculate_forward(size_t N, size_t H, size_t ndof);
    static Matrix_t calculate_backward(size_t N, size_t H, size_t ndof);

    DFT() = default;
    DFT(size_t H, size_t N, size_t ndof = 1);    

};

} // namespace math