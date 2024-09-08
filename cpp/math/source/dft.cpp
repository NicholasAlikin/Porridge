#include "dft.hpp"

namespace math {


Matrix_t DFT::calculate_forward_basic(size_t N, size_t H) {
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
Matrix_t DFT::calculate_backward_basic(size_t N, size_t H) {
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

Matrix_t DFT::calculate_forward(size_t N, size_t H, size_t ndof) {
    return kron(DFT::calculate_forward_basic(N,H), eye<double>(ndof));
}
Matrix_t DFT::calculate_backward(size_t N, size_t H, size_t ndof) {
    return kron(DFT::calculate_backward_basic(N,H), eye<double>(ndof));
}

DFT::DFT(size_t H, size_t N, size_t ndof)
        : H(H), N(N), ndof(ndof) {

}


} // namespace math