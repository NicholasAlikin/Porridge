#include "dft.hpp"

namespace npath {

size_t DFT::frequency_basic_size() const {
    return 2*H+1;
}
size_t DFT::frequency_size() const {
    return frequency_basic_size()*ndof;
}

size_t DFT::time_basic_size() const {
    return N;
}
size_t DFT::time_size() const {
    return time_basic_size()*ndof;
}

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

math::vector_t<double,3> DFT::calculate_forward_matrix_basic(size_t H, size_t N) {
    math::vector_t<double,3> M_forward(N,math::vector_t<double,2>(2*H+1));
    math::vector<double> l(2*H+1);
    l[0] = 1.0;

    auto itl_end = l.end();
    decltype(itl_end) itl;
    decltype(M_forward[0].begin()) itM;

    for (size_t n = 0; n < N; ++n) {
        // store l vector
        for (size_t h = 1; h <= H; ++h) {
            l[2*h-1] = std::cos(h* 2*math::pi/N * n);
            l[2*h]   = std::sin(h* 2*math::pi/N * n);
        }
        M_forward[n][0] = 0.5*l;
        itl = l.begin() + 1;
        itM = M_forward[n].begin()+1;
        while (itl < itl_end) {
            *itM = (*itl)*l;
            ++itl;
            ++itM;
        }
    }

    return M_forward*(2./N);
}

DFT::DFT(size_t H, size_t N, size_t ndof)
        : H(H), N(N), ndof(ndof)
        , forward_basic(DFT::calculate_forward_basic(H,N))
        , backward_basic(DFT::calculate_backward_basic(H,N))
        , derivative_basic(DFT::calculate_derivative_basic(H,N))
        , derivative2_basic(DFT::calculate_derivative2_basic(derivative_basic))
{}


DFT::DFT(const DFT& other)
        : H(other.H), N(other.N), ndof(other.ndof)
        , forward_basic(other.forward_basic)
        , backward_basic(other.backward_basic)
        , derivative_basic(other.derivative_basic)
        , derivative2_basic(other.derivative2_basic)

        , forward_matrix(other.forward_matrix)
        , forward_matrix_dot_derivative(other.forward_matrix_dot_derivative)
        , forward_matrix_dot_derivative2(other.forward_matrix_dot_derivative2)

        , mass_freq(other.mass_freq)
        , damp_freq(other.damp_freq)
        , stif_freq(other.stif_freq)

        , derivative_forward_helper(other.derivative_forward_helper)
        , lb0_helper(other.lb0_helper)
{}


DFT::DFT(DFT&& other)
        : H(other.H), N(other.N), ndof(other.ndof)
        , forward_basic(std::move(other.forward_basic))
        , backward_basic(std::move(other.backward_basic))
        , derivative_basic(std::move(other.derivative_basic))
        , derivative2_basic(std::move(other.derivative2_basic))

        , forward_matrix(std::move(other.forward_matrix))
        , forward_matrix_dot_derivative(std::move(other.forward_matrix_dot_derivative))
        , forward_matrix_dot_derivative2(std::move(other.forward_matrix_dot_derivative2))

        , mass_freq(std::move(other.mass_freq))
        , damp_freq(std::move(other.damp_freq))
        , stif_freq(std::move(other.stif_freq))

        , derivative_forward_helper(std::move(other.derivative_forward_helper))
        , lb0_helper(std::move(other.lb0_helper))
{}


math::vector<double> DFT::derivative(math::vector<double>&& vec) const {
    return derivative(std::move(vec), ndof);
}


math::vector<double> DFT::derivative(math::vector<double>&& vec, size_t ndofs) const {
    math::vector<double> res = std::move(vec);
    __derivative(res, ndofs);
    return res;
}



math::vector<double> DFT::derivative2(math::vector<double>&& vec) const {
    math::vector<double> res = std::move(vec);
    __derivative2(res);
    return res;
}


void DFT::jac_time_domain_precomp(size_t base_ndof) {
    forward_matrix = DFT::calculate_forward_matrix_basic(H,N);
    forward_matrix_dot_derivative  = math::zeros<double>(forward_matrix);
    forward_matrix_dot_derivative2 = math::zeros<double>(forward_matrix);

    for (size_t t = 0; t < N; ++t) {
        forward_matrix_dot_derivative[t] = math::dot(forward_matrix[t],
                                                     derivative_basic);

        forward_matrix_dot_derivative2[t] = math::dot(forward_matrix[t],
                                                      derivative2_basic);
    }
    
    mass_freq = math::zeros<double>(base_ndof*frequency_basic_size(),
                                    base_ndof*frequency_basic_size());
    damp_freq = math::zeros<double>(base_ndof*frequency_basic_size(),
                                    base_ndof*frequency_basic_size());
    stif_freq = math::zeros<double>(base_ndof*frequency_basic_size(),
                                    base_ndof*frequency_basic_size());
    
    /**/
    derivative_forward_helper = math::zeros<double>(N,frequency_basic_size());
    math::vector<double> lb(2*H+1);
    lb[0] = 0.5;

    auto itlb_end = lb.end();
    decltype(itlb_end) itl;
    
    for (size_t n = 0; n < N; ++n) { 
        // store l vector
        for (size_t h = 1; h <= H; ++h) {
            lb[2*h-1] = std::cos(h* 2.*math::pi/N * n);
            lb[2*h]   = std::sin(h* 2.*math::pi/N * n);
        }
        if (n == 0) {
            lb0_helper = lb*(-2.);
            continue; // derivative_forward_helper[n = 0] == 0
        }
        math::dot(derivative_basic,lb*(2./N* 2.*math::pi*n/N)
                  ,derivative_forward_helper[n]);
    }
}

void DFT::jac_time_domain_inc(math::vector<math::vector_slice<double>>& matrix
                            ,const math::vector_t<double,2>& mass
                            ,const math::vector_t<double,2>& damp
                            ,const math::vector_t<double,2>& stif
                            ,double freq, size_t time) const
{
    if (time == 0) {
        math::fill(mass_freq,0.0);
        math::fill(damp_freq,0.0);
        math::fill(stif_freq,0.0);
    }

    math::addkron(forward_matrix_dot_derivative2[time], mass, mass_freq);
    math::addkron(forward_matrix_dot_derivative[time], damp, damp_freq);
    math::addkron(forward_matrix[time], stif, stif_freq);
    

    /* last loop iteration */
    if (time == N-1) {
        matrix = (mass_freq*(freq*freq) + damp_freq*(freq) + stif_freq);
    }
}

    
void DFT::jac_extended_time_domain_inc(math::vector<math::vector_slice<double>>& drdq
                                ,const math::vector_t<double,2>& mass
                                ,const math::vector_t<double,2>& damp
                                ,const math::vector_t<double,2>& stif
                                ,double freq, size_t time
                                     , math::vector_slice<double>& drdw
                                ,const math::vector_const_slice<double>& q
                                ,const math::vector_const_slice<double>& r) const
{
    if (time == 0) {
        math::fill(mass_freq,0.0);
        math::fill(damp_freq,0.0);
        math::fill(stif_freq,0.0);
    }

    math::addkron(forward_matrix_dot_derivative2[time], mass, mass_freq);
    math::addkron(forward_matrix_dot_derivative[time], damp, damp_freq);
    math::addkron(forward_matrix[time], stif, stif_freq);
    

    /* last loop iteration */
    if (time == N-1) {
        drdq = (mass_freq*(freq*freq) + damp_freq*(freq) + stif_freq);
        math::dot((mass_freq*(freq*2) + damp_freq), q, drdw);
    }

}


void DFT::time_domain(const math::vector<double>& q, double freq
                          ,math::vector<double>& u
                          ,math::vector<double>& dudt
                          ,math::vector<double>& d2udt2, size_t ndofs) const
{   
    math::fill(u,0.0);
    math::fill(dudt,0.0);
    math::fill(d2udt2,0.0);
    backward(q, u, ndofs);
    
    math::vector<double> subq = q*freq;
    subq = derivative(std::move(subq), ndofs);
    backward(subq, dudt, ndofs);
    
    subq *= freq;
    subq = derivative(std::move(subq), ndofs);
    backward(subq, d2udt2, ndofs);
}


} // namespace npath