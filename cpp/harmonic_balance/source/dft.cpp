#include "dft.hpp"

namespace npath {


math::vector_t<double,2> DFT::calculate_forward_basic(size_t H, size_t N) {
    math::vector_t<double,2> T_forward = math::zeros<double>(2*H+1, N);
    T_forward[0] = math::repmat<double>(1.0/N, N);
    
    size_t h, n;
    for (h = 1; h <= H; ++h) {
        for (n = 0; n < N; ++n) {
            T_forward[2*h-1][n] = 2.0/N* std::cos(h* 2*math::pi/N * n);
            T_forward[2*h][n]   = 2.0/N* std::sin(h* 2*math::pi/N * n);
        }
    }
    return T_forward;
}
math::vector_t<double,2> DFT::calculate_backward_basic(size_t H, size_t N) {
    math::vector_t<double,2> T_backward = math::zeros<double>(N, 2*H+1);
    
    size_t h, n;
    for (n = 0; n < N; ++n) {
        T_backward[n][0] = 1.0;
        for (h = 1; h <= H; ++h) {
            T_backward[n][2*h-1] = std::cos(h* 2*math::pi/N * n);
            T_backward[n][2*h]   = std::sin(h* 2*math::pi/N * n);
        }
    }
    return T_backward;
}

math::vector_t<double,2> DFT::calculate_forward(size_t H, size_t N, size_t ndof) {
    return math::kron(DFT::calculate_forward_basic(H,N), math::eye<double>(ndof));
}
math::vector_t<double,2> DFT::calculate_backward(size_t H, size_t N, size_t ndof) {
    return math::kron(DFT::calculate_backward_basic(H,N), math::eye<double>(ndof));
}

math::vector_t<double,2> DFT::calculate_derivative_basic(size_t H, size_t) {
    math::vector_t<double,2> T_derivative = math::zeros<double>(2*H+1, 2*H+1);
    double k = 1;
    for(size_t h = 1; h < H+1; ++h, ++k) {
        T_derivative[2*h-1][2*k  ] =  k;
        T_derivative[2*h  ][2*k-1] = -k;
    }
    return T_derivative;
}

math::vector_t<double,2> DFT::calculate_derivative2_basic(size_t H, size_t N) {
    math::vector_t<double,2> T_derivative2 = DFT::calculate_derivative_basic(H, N);
    return dot(T_derivative2,T_derivative2);
}

math::vector_t<double,2> DFT::calculate_derivative2_basic(const math::vector_t<double,2>& T_derivative2_basic) {
    return dot(T_derivative2_basic,T_derivative2_basic);
}

math::vector_t<double,2> DFT::calculate_derivative(size_t H, size_t N, size_t ndof) {
    math::vector_t<double,2> T_derivative = math::kron(DFT::calculate_derivative_basic(H,N), math::eye<double>(ndof));
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


math::vector<double> DFT::transform(const math::vector_t<double,2>& transform_matrix, size_t res_size, const math::vector<double>& vec) const {
    
    math::vector<double> res = math::zeros<double>(res_size);
    math::Slice<decltype(res.begin()),decltype(res.end())> sl_res;
    math::Slice<decltype(vec.begin()),decltype(vec.end())> sl_vec;
    
    for (size_t dof = 0; dof < ndof; ++dof){
        sl_res.new_slice(res.begin()+dof, res.end(), ndof);
        sl_vec.new_slice(vec.begin()+dof, vec.end(), ndof);
        math::dot(transform_matrix,sl_vec,sl_res);
    }
    return res;
}

math::vector<double> DFT::forward(const math::vector<double>& vec) const {
    /*vec.size() == N ndof, res.size() == (2H+1)ndof*/
    return transform(forward_basic,(2*H+1)*ndof,vec);
}

math::vector<double> DFT::backward(const math::vector<double>& vec) const {
    /*vec.size() == (2H+1)ndof, res.size() == N ndof*/
    return transform(backward_basic,N*ndof,vec);
}

void DFT::__derivative(math::vector<double>& res) const {
    auto res_dof = res.begin();
    math::Slice<decltype(res.begin()),decltype(res.end())> sl_cos, sl_sin;
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


math::vector<double> DFT::derivative(const math::vector<double>& vec) const {
    math::vector<double> res = vec;
    __derivative(res);
    return res;
}

math::vector<double> DFT::derivative(math::vector<double>&& vec) const {
    math::vector<double> res = std::move(vec);
    __derivative(res);
    return res;
}

void DFT::__derivative2(math::vector<double>& res) const {
    auto res_dof = res.begin();
    math::Slice<decltype(res.begin()),decltype(res.end())> sl_cos, sl_sin;
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
math::vector<double> DFT::derivative2(const math::vector<double>& vec) const {
    math::vector<double> res = vec;
    __derivative2(res);
    return res;
}

math::vector<double> DFT::derivative2(math::vector<double>&& vec) const {
    math::vector<double> res = std::move(vec);
    __derivative2(res);
    return res;
}
} // namespace npath