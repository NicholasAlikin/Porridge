#pragma once
#include "linalg.hpp"

namespace npath {

class DFT {

public:
    size_t H = 0;
    size_t N = 0;
    size_t ndof = 1;

    // vector transform for ndof = 1
    math::vector_t<double,2> forward_basic;
    math::vector_t<double,2> backward_basic;
    math::vector_t<double,2> derivative_basic;  // for ndof = 1
    math::vector_t<double,2> derivative2_basic; // for ndof = 1
    // matrix transform for ndof = 1
// private:
    mutable math::vector_t<double,3> forward_matrix;
    mutable math::vector_t<double,3> forward_matrix_dot_derivative;
    mutable math::vector_t<double,3> forward_matrix_dot_derivative2;

    /* Tangent matrices in frequency domain */
    mutable math::vector_t<double,2> mass_freq;
    mutable math::vector_t<double,2> damp_freq;
    mutable math::vector_t<double,2> stif_freq;

    /* Other */
    mutable math::vector_t<double,2> derivative_forward_helper;
    mutable math::vector<double> lb0_helper;


public:
    

    size_t frequency_size() const;
    size_t frequency_basic_size() const;
    size_t time_size() const;
    size_t time_basic_size() const;

    static math::vector_t<double,2> calculate_forward_basic(size_t H, size_t N);
    static math::vector_t<double,2> calculate_backward_basic(size_t H, size_t N);
    static math::vector_t<double,2> calculate_derivative_basic(size_t H, size_t N);
    static math::vector_t<double,2> calculate_derivative2_basic(size_t H, size_t N);
    static math::vector_t<double,2> calculate_derivative2_basic(const math::vector_t<double,2>& T_derivative_basic);

    static math::vector_t<double,3> calculate_forward_matrix_basic(size_t H, size_t N);

    DFT() = default;
    DFT(size_t H, size_t N, size_t ndof = 1);
    DFT(const DFT& other);
    DFT(DFT&& other);

private:
    template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
    void transform(const math::vector_t<double,2>& transform_matrix,
                   const V1& vec, V2& res, size_t ndofs) const;
    

    template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
    void transform(const math::vector_t<double,2>& transform_matrix,
                   const V1& vec, V2& res) const;


    template <math::Vector V>
    void __derivative(V& res, size_t ndofs) const;
    

    template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
    void __derivative(const V1& vec, V2& res, size_t ndofs) const;
    
    
    template <math::Vector V>
    void __derivative2(V& res) const;
    

    template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
    void __derivative2(const V1& vec, V2& res) const;

public:
    
    template <math::Vector V>
    math::vector<double> forward(const V& vec) const;
    
    
    template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
    void                 forward(const V1& vec, V2& res, size_t ndofs) const;


    template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
    void                 forward(const V1& vec, V2& res) const;
    

    template <math::Vector V>
    math::vector<double> backward(const V& vec) const;


    template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
    void                 backward(const V1& vec, V2& res, size_t ndofs) const;


    template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
    void                 backward(const V1& vec, V2& res) const;


    template <math::Vector V>
    math::vector<double> derivative(const V& vec) const;

    template <math::Vector V>
    math::vector<double> derivative(const V& vec, size_t ndofs) const;


    template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
    void                 derivative(const V1& vec, V2& res) const;


    template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
    void                 derivative(const V1& vec, V2& res, size_t ndofs) const;


    math::vector<double> derivative(math::vector<double>&& vec) const;


    math::vector<double> derivative(math::vector<double>&& vec, size_t ndofs) const;
    

    template <math::Vector V>
    math::vector<double> derivative2(const V& vec) const;


    template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
    void                 derivative2(const V1& vec, V2& res) const;

    math::vector<double> derivative2(math::vector<double>&& vec) const;
    

    void jac_time_domain_precomp(size_t base_ndof);
    void jac_time_domain_inc(math::vector<math::vector_slice<double>>& matrix
                            ,const math::vector_t<double,2>& mass
                            ,const math::vector_t<double,2>& damp
                            ,const math::vector_t<double,2>& stif
                            , double freq, size_t time) const;

    /* Iteration of calculation Jacobi matrix of system in frequency domain as
        1. Transfrorm state vector from frequency to time domain,
        2. Calculate tangent matrices `mass`, `damp` and `stif` in time domain
        3. For each `time` sample call this method for the next iteration
    Also iteration of calculation ... */
    void jac_extended_time_domain_inc(math::vector<math::vector_slice<double>>& drdq
                               ,const math::vector_t<double,2>& mass
                               ,const math::vector_t<double,2>& damp
                               ,const math::vector_t<double,2>& stif
                               ,double freq, size_t time
                                    , math::vector_slice<double>& drdw
                               ,const math::vector_const_slice<double>& q
                               ,const math::vector_const_slice<double>& r) const;

    void time_domain(const math::vector<double>& q, double freq
                          ,math::vector<double>& u
                          ,math::vector<double>& dudt
                          ,math::vector<double>& d2udt2, size_t ndof) const;


};



template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
void DFT::transform(const math::vector_t<double,2>& transform_matrix
                  , const V1& vec, V2& res, size_t ndofs) const
{
    math::Slice sl_res(res.begin(), res.end(), ndofs);
    math::Slice sl_vec(vec.begin(), vec.end(), ndofs);
    
    for (size_t dof = 0; dof < ndofs; ++dof) {
        sl_res.update_from(res.begin()+dof); // slices size remains constant
        sl_vec.update_from(vec.begin()+dof);
        math::dot(transform_matrix,sl_vec,sl_res);
    }
}


template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
void DFT::transform(const math::vector_t<double,2>& transform_matrix
                  , const V1& vec, V2& res) const
{
    transform(transform_matrix, vec, res, ndof);
}


template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
void DFT::forward(const V1& vec, V2& res, size_t ndofs) const {
    transform(forward_basic, vec, res, ndofs);
}


template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
void DFT::forward(const V1& vec, V2& res) const {
    forward(vec, res, ndof);
}


template <math::Vector V>
math::vector<double> DFT::forward(const V& vec) const {
    auto res = math::zeros<double>(frequency_size());
    forward(vec,res);
    return res;
}


template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
void DFT::backward(const V1& vec, V2& res, size_t ndofs) const {
    transform(backward_basic,vec,res,ndofs);
}


template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
void DFT::backward(const V1& vec, V2& res) const {
    backward(vec,res,ndof);
}


template <math::Vector V>
math::vector<double> DFT::backward(const V& vec) const {
    auto res = math::zeros<double>(time_size());
    backward(vec,res);
    return res;
}



template <math::Vector V>
void DFT::__derivative(V& res, size_t ndofs) const {
    auto res_dof = res.begin();
    math::Slice sl_cos(res_dof+1*ndofs, res.end(), 2*ndofs)
              , sl_sin(res_dof+2*ndofs, res.end(), 2*ndofs);
    decltype(sl_cos.begin()) it_cos,it_sin;
    size_t h;
    double tmp;
    for (size_t dof = 0; dof < ndofs; ++dof) {
        *res_dof = 0;
        h = 1;
        sl_cos.update_from(res_dof+1*ndofs);
        sl_sin.update_from(res_dof+2*ndofs);
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


template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
void DFT::__derivative(const V1& vec, V2& res, size_t ndofs) const {
    auto vec_dof = vec.begin();
    auto res_dof = res.begin();

    math::Slice sl0_cos(vec_dof+1*ndofs, vec.end(), 2*ndofs);
    math::Slice sl_cos( res_dof+1*ndofs, res.end(), 2*ndofs)
              , sl_sin( res_dof+2*ndofs, res.end(), 2*ndofs);
    
    decltype(sl0_cos.begin()) it0_cos,it0_sin;
    decltype(sl_cos.begin()) it_cos,it_sin;
    
    size_t h;
    for (size_t dof = 0; dof < ndofs; ++dof) {
        *res_dof = 0;
        
        sl0_cos.update_from(vec_dof+1*ndofs);
        it0_cos = sl0_cos.begin();
        
        sl_cos.update_from(res_dof+1*ndofs);
        sl_sin.update_from(res_dof+2*ndofs);
        it_cos = sl_cos.begin();
        it_sin = sl_sin.begin();

        for (h = 1; h <= H; ++h
                            ,++it0_cos
                            ,++it_cos
                            ,++it_sin)
        {
            *it_cos =  (*it_sin )*h;
            *it_sin = -(*it0_cos)*h;
        }
        ++res_dof;
        ++vec_dof;
    }
}


template <math::Vector V>
math::vector<double> DFT::derivative(const V& vec) const {
    math::vector<double> res = vec;
    __derivative(vec,res,ndof);
    return res;
}


template <math::Vector V>
math::vector<double> DFT::derivative(const V& vec, size_t ndofs) const {
    math::vector<double> res = vec;
    __derivative(vec,res,ndofs);
    return res;
}


template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
void DFT::derivative(const V1& vec, V2& res) const {
    derivative(vec,res,ndof);
}


template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
void DFT::derivative(const V1& vec, V2& res, size_t ndofs) const {
    __derivative(vec,res,ndofs);
}


template <math::Vector V>
void DFT::__derivative2(V& res) const {
    auto res_dof = res.begin();
    math::Slice sl_cos(res_dof+1*ndof, res.end(), 2*ndof)
              , sl_sin(res_dof+2*ndof, res.end(), 2*ndof);
    decltype(sl_cos.begin()) it_cos,it_sin;
    size_t h;
    for (size_t dof = 0; dof < ndof; ++dof) {
        *res_dof = 0;
        h = 1;
        sl_cos.update_from(res_dof+1*ndof);
        sl_sin.update_from(res_dof+2*ndof);
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



template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
void DFT::__derivative2(const V1& vec, V2& res) const {
    auto vec_dof = vec.begin();
    auto res_dof = res.begin();
    
    math::Slice sl0_cos(vec_dof+1*ndof, vec.end(), 2*ndof)
              , sl0_sin(vec_dof+2*ndof, vec.end(), 2*ndof);
    math::Slice sl_cos( res_dof+1*ndof, res.end(), 2*ndof)
              , sl_sin( res_dof+2*ndof, res.end(), 2*ndof);
    
    decltype(sl0_cos.begin()) it0_cos,it0_sin;
    decltype(sl_cos.begin() ) it_cos, it_sin;
    size_t h;
    for (size_t dof = 0; dof < ndof; ++dof) {
        *res_dof = 0;
        
        sl0_cos.update_from(vec_dof+1*ndof);
        sl0_sin.update_from(vec_dof+2*ndof);
        it0_cos = sl0_cos.begin();
        it0_sin = sl0_sin.begin();

        sl_cos.update_from(res_dof+1*ndof);
        sl_sin.update_from(res_dof+2*ndof);
        it_cos = sl_cos.begin();
        it_sin = sl_sin.begin();
        
        for (h = 1; h <= H; ++h
                           ,++it0_cos
                           ,++it0_sin
                           ,++it_cos
                           ,++it_sin)
        {
            *it_cos = -h*h* (*it0_cos);
            *it_sin = -h*h* (*it0_sin);
        }

        ++vec_dof;
        ++res_dof;
    }
}


template <math::Vector V>
math::vector<double> DFT::derivative2(const V& vec) const {
    math::vector<double> res = vec;
    __derivative2(vec,res);
    return res;
}


template <math::Vector V1, math::ArithmeticVectorsLike<V1> V2>
void DFT::derivative2(const V1& vec, V2& res) const {
    __derivative2(vec,res);
}

} // namespace npath