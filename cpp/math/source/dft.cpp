#include "dft.hpp"

namespace math {


Matrix_t DFT::calculate_forward_basic(size_t H, size_t N) {
    Matrix_t T_forward = zeros<double>(2*H+1, N);
    T_forward[0] = repmat<double>(1.0/N, N);
    
    size_t h, n;
    for (h = 1; h <= H; ++h) {
        for (n = 0; n < N; ++n) {
            T_forward[2*h-1][n] = 2.0/N* std::cos(h* 2*pi/N * n);
            T_forward[2*h][n]   = 2.0/N* std::sin(h* 2*pi/N * n);
        }
    }
    return T_forward;
}
Matrix_t DFT::calculate_backward_basic(size_t H, size_t N) {
    Matrix_t T_backward = zeros<double>(N, 2*H+1);
    
    size_t h, n;
    for (n = 0; n < N; ++n) {
        T_backward[n][0] = 1.0;
        for (h = 1; h <= H; ++h) {
            T_backward[n][2*h-1] = std::cos(h* 2*pi/N * n);
            T_backward[n][2*h]   = std::sin(h* 2*pi/N * n);
        }
    }
    return T_backward;
}

Matrix_t DFT::calculate_forward(size_t H, size_t N, size_t ndof) {
    return kron(DFT::calculate_forward_basic(H,N), eye<double>(ndof));
}
Matrix_t DFT::calculate_backward(size_t H, size_t N, size_t ndof) {
    return kron(DFT::calculate_backward_basic(H,N), eye<double>(ndof));
}

Matrix_t DFT::calculate_derivative_basic(size_t H, size_t N) {
    Matrix_t T_derivative = zeros<double>(2*H+1, 2*H+1);
    double k = 1;
    for(size_t h = 1; h < H+1; ++h, ++k) {
        T_derivative[2*h-1][2*k  ] =  k;
        T_derivative[2*h  ][2*k-1] = -k;
    }
    return T_derivative;
}

Matrix_t DFT::calculate_derivative2_basic(size_t H, size_t N) {
    Matrix_t T_derivative2 = DFT::calculate_derivative_basic(N, H);
    return dot(T_derivative2,T_derivative2);
}

Matrix_t DFT::calculate_derivative2_basic(const Matrix_t& T_derivative2_basic) {
    return dot(T_derivative2_basic,T_derivative2_basic);
}

Matrix_t DFT::calculate_derivative(size_t H, size_t N, size_t ndof) {
    Matrix_t T_derivative = kron(DFT::calculate_derivative_basic(H,N), eye<double>(ndof));
    return T_derivative;
}

DFT::DFT(size_t H, size_t N, size_t ndof)
        : H(H), N(N), ndof(ndof)
        // , forward_basic(DFT::calculate_forward_basic(H,N))
        , derivative_basic(DFT::calculate_derivative_basic(H,N))
        , derivative2_basic(DFT::calculate_derivative2_basic(derivative_basic))
        , forward(DFT::calculate_forward(H,N,ndof))
        , backward(DFT::calculate_backward(H,N,ndof))
        , derivative(DFT::calculate_derivative(H,N,ndof))
        // , derivative2(DFT::calculate_derivative2(H,N,ndof))
        {
}


} // namespace math