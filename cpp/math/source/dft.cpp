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
        , forward_basic(DFT::calculate_forward_basic(H,N))
        , backward_basic(DFT::calculate_backward_basic(H,N))
        , derivative_basic(DFT::calculate_derivative_basic(H,N))
        , derivative2_basic(DFT::calculate_derivative2_basic(derivative_basic))
        {
}


Vector_t DFT::transform(const Matrix_t& transform_matrix, size_t res_size, const Vector_t& vec) const {
    
    Vector_t res = zeros<double>(res_size);
    Slice<decltype(res.begin()),decltype(res.end())> sl_res;
    Slice<decltype(vec.begin()),decltype(vec.end())> sl_vec;
    
    for (size_t dof = 0; dof < ndof; ++dof){
        sl_res.new_slice(res.begin()+dof, res.end(), ndof);
        sl_vec.new_slice(vec.begin()+dof, vec.end(), ndof);
        dot(transform_matrix,sl_vec,sl_res);
    }
    return res;
}

Vector_t DFT::forward(const Vector_t& vec) const {
    /*vec.size() == N ndof, res.size() == (2H+1)ndof*/
    return transform(forward_basic,(2*H+1)*ndof,vec);
}

Vector_t DFT::backward(const Vector_t& vec) const {
    /*vec.size() == (2H+1)ndof, res.size() == N ndof*/
    return transform(backward_basic,N*ndof,vec);
}

void DFT::__derivative(Vector_t& res) const {
    auto res_dof = res.begin();
    Slice<decltype(res.begin()),decltype(res.end())> sl_cos, sl_sin;
    decltype(sl_cos.begin()) it_cos,it_sin;
    size_t h;
    double tmp;
    for (size_t dof = 0; dof < ndof; ++dof) {
        *res_dof = 0;
        h = 1;
        sl_cos.new_slice(res_dof+1*ndof, res.end(), 2*ndof);
        sl_sin.new_slice(res_dof+2*ndof, res.end(), 2*ndof);
        it_cos = sl_cos.begin();
        it_sin = sl_sin.begin();
        while (h <= H) {
            tmp = *it_cos;
            *it_cos = (*it_sin)*h;
            *it_sin = -tmp*h;
            ++h; ++it_cos; ++it_sin;
        }
        ++res_dof;
    }
}


Vector_t DFT::derivative(const Vector_t& vec) const {
    Vector_t res = vec;
    __derivative(res);
    return res;
}

Vector_t DFT::derivative(Vector_t&& vec) const {
    Vector_t res = std::move(vec);
    __derivative(res);
    return res;
}

void DFT::__derivative2(Vector_t& res) const {
    auto res_dof = res.begin();
    Slice<decltype(res.begin()),decltype(res.end())> sl_cos, sl_sin;
    decltype(sl_cos.begin()) it_cos,it_sin;
    size_t h;
    for (size_t dof = 0; dof < ndof; ++dof) {
        *res_dof = 0;
        h = 1;
        sl_cos.new_slice(res_dof+1*ndof, res.end(), 2*ndof);
        sl_sin.new_slice(res_dof+2*ndof, res.end(), 2*ndof);
        it_cos = sl_cos.begin();
        it_sin = sl_sin.begin();
        while (h <= H) {
            (*it_cos) *= -h*h;
            (*it_sin) *= -h*h;
            ++h; ++it_cos; ++it_sin;
        }
        ++res_dof;
    }
}
Vector_t DFT::derivative2(const Vector_t& vec) const {
    Vector_t res = vec;
    __derivative2(res);
    return res;
}

Vector_t DFT::derivative2(Vector_t&& vec) const {
    Vector_t res = std::move(vec);
    __derivative2(res);
    return res;
}
} // namespace math