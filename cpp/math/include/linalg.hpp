#pragma once

#include "vector.hpp"
#include "slice.hpp"
#include "vectorlike_functions.hpp"
#include "vectorlike_type_traits.hpp"
#include <cmath>
#include <numbers>
#include <string>

namespace math {

constexpr double pi = std::numbers::pi;

/*==================================
// Generate vector-like objects
==================================*/
namespace detail {
// zeros - vector-like full of zeros (default constructed vector::basic_value_type)
template <typename T, size_t Dim>
requires (Dim == 1)
vector_t<T,Dim> zeros_helper(size_t sz) {
    return vector_t<T,Dim>(sz);
}

template <typename T, size_t Dim>
vector_t<T,Dim> zeros_helper(size_t sz, auto... szs) {
    return vector_t<T,Dim>(sz, zeros_helper<T,Dim-1>(szs...));
}

template <typename T, size_t Dim>
requires (Dim == 1)
vector_t<T,Dim> zeros_helper(const size_t* sizes) {
    return vector_t<T,Dim>(*sizes);
}

template <typename T, size_t Dim>
vector_t<T,Dim> zeros_helper(const size_t* sizes) {
    size_t sz = *sizes;
    return vector_t<T,Dim>(sz, zeros_helper<T,Dim-1>(++sizes));
}
} // namespace detail

template <typename T, typename... Sizes>
requires (std::is_convertible_v<Sizes,size_t> && ...)
auto zeros(Sizes... sizes) 
        -> vector_t<T,sizeof...(Sizes)> {
    return detail::zeros_helper<T,sizeof...(Sizes)>(sizes...);
}

template <typename T, size_t... sizes>
auto zeros(std::integer_sequence<size_t,sizes...> int_seq) 
        -> vector_t<T,sizeof...(sizes)> {
    return detail::zeros_helper<T,int_seq.size()>((sizes, int_seq.size()) ...);
}

template <typename T, size_t Dim>
vector_t<T,Dim> zeros(const std::array<size_t,Dim>& sizes) {
    return detail::zeros_helper<T,Dim>(sizes.begin());
}

template <typename T, VectorLike Vec>
auto zeros(const Vec& v)
        -> vector_t<T,vector_dim_v<Vec>> {
    return zeros<T>(size(v));
}

// repmat - vector-like full of given value with type vector::basic_value_type
namespace detail {
template <typename T, size_t Dim>
requires (Dim == 1)
vector_t<T,Dim> repmat_helper(const T& value, size_t sz) {
    return vector_t<T,Dim>(sz,value);
}
template <typename T, size_t Dim>
vector_t<T,Dim> repmat_helper(const T& value, size_t sz, auto... sizes) {
    return vector_t<T,Dim>(sz, repmat_helper<T,Dim-1>(value, sizes...));
}

template <typename T, size_t Dim>
requires (Dim == 1)
vector_t<T,Dim> repmat_helper(const T& value, const size_t* sizes) {
    return vector_t<T,Dim>(*sizes, value);
}

template <typename T, size_t Dim>
vector_t<T,Dim> repmat_helper(const T& value, const size_t* sizes) {
    const size_t sz = *sizes;
    return vector_t<T,Dim>(sz, detail::repmat_helper<T,Dim-1>(value, ++sizes));
}
} // namespace detail

template <typename T, typename... Sizes>
requires (std::is_convertible_v<Sizes,size_t> && ...)
vector_t<T,sizeof...(Sizes)> repmat(const T& value, Sizes... sizes) {
    return detail::repmat_helper<T,sizeof...(Sizes)>(value, sizes...);
}

template <typename T, size_t Dim>
vector_t<T,Dim> repmat(const T& value, std::array<size_t,Dim> sizes) {
    return detail::repmat_helper<T,Dim>(value, sizes.begin());
}

template <typename T, VectorLike Vec>
auto repmat(const T& value, const Vec& v)
        -> vector_t<T,vector_dim_v<Vec>> {
    return repmat<T>(value, size(v));
}

// eye - vector_t<T,2> equal to zeros except of diagonal elements
// which equal to 1
template <typename T>
vector_t<T,2> eye(size_t sz) {
    vector_t<T,2> v = zeros<T>(sz,sz);
    for (size_t i = 0; i < sz; ++i) {
        v[i][i] = 1;
    }
    return v;
}


template <Matrix M>
auto det_(const M& vec)
        -> typename M::basic_value_type {
    const size_t sz = vec.size();
    
    if (vec.size() == 1)
        return vec[0][0];
    if (vec.size() == 2)
        return vec[0][0]*vec[1][1] - vec[0][1]*vec[1][0];

    int j = 1;
    size_t k,l,m,n;
    typename M::basic_value_type res = 0;
    auto temp = zeros<typename M::basic_value_type>(sz-1,sz-1);
    for (size_t i = 0; i < sz; ++i, j*=-1) {
        
        for (k=0,m=1; m < sz; ++m) {
            for (l=0,n=0; n < sz; ++n) {
                if (n == i) continue;
                temp[k][l] = vec[m][n];
                ++l;
            }
            ++k;
        }

        res += vec[0][i]*j*det_(temp);
    }
    return res;

}

// double det_(const vector_t<double,2>& mat);

template <Matrix M>
auto det(const M& mat)
        -> typename M::basic_value_type {
    auto sz = size(mat);
    if (sz[0] != sz[1]) {
        throw std::logic_error("Cannot calculate matrix determinant! Matrix dimensions must be equal!");
    }

    return det_(mat);
}

template <Vector V>
auto sum(const V& vec)
        -> typename V::basic_value_type
{
    typename V::basic_value_type res = 0;
    auto it = vec.begin(), end = vec.end();
    while (it < end) {
        res += *it;
        ++it;
    }
    return res;
}

template <Vector V>
auto product(const V& vec)
        -> typename V::basic_value_type
{
    typename V::basic_value_type res = 1;
    auto it = vec.begin(), end = vec.end();
    while (it < end) {
        res *= *it;
        ++it;
    }
    return res;
}

namespace detail {

template <Matrix M1, ArithmeticVectorsLike<M1> M2>
void check_dot_operands_sizes(const M1& mat1, const M2& mat2) {
    auto sz1 = size(mat1);
    auto sz2 = size(mat2);
    if (sz1[1] != sz2[0]) {
        throw std::logic_error(std::string("Cannot calculate dot(mat[")
                             + std::to_string(sz1[0])
                             + std::string(",")
                             + std::to_string(sz1[1])
                             + std::string("], mat[")
                             + std::to_string(sz2[0])
                             + std::string(",")
                             + std::to_string(sz2[1])
                             + std::string("])! Incorrect matrices shapes!"));
    }
}

template <VectorLike V3D, Vector V>
requires Matrix<typename V3D::value_type>
void check_dot_operands_sizes(const V3D& v3d, const V& v) {
    auto sz1 = size(v3d);
    auto sz2 = size(v);
    if (sz1[2] != sz2[0]) {
        throw std::logic_error(std::string("Cannot calculate dot(vec3d[")
                             + std::to_string(sz1[0])
                             + std::string(",")
                             + std::to_string(sz1[1])
                             + std::string(",")
                             + std::to_string(sz1[2])
                             + std::string("], vec[")
                             + std::to_string(sz2[0])
                             + std::string("])! Incorrect 3d or 1d vectors shapes!"));
    }
}

template <VectorLike V3D, Vector V>
requires Matrix<typename V3D::value_type>
void check_dot_operands_sizes(const V& v, const V3D& v3d) {
    auto sz1 = size(v);
    auto sz2 = size(v3d);
    if (sz1[0] != sz2[0]) {
        throw std::logic_error(std::string("Cannot calculate dot(vec[")
                             + std::to_string(sz1[0])
                             + std::string("], vec3d[")
                             + std::to_string(sz2[0])
                             + std::string(",")
                             + std::to_string(sz2[1])
                             + std::string(",")
                             + std::to_string(sz2[2])
                             + std::string("])! Incorrect 1d or 3d vectors shapes!"));
    }
}

template <Matrix M, Vector V>
requires HaveGeneralType<M,V>
void check_dot_operands_sizes(const M& mat, const V& vec) {
    auto sz = size(mat);
    if (sz[1] != vec.size()) {
        throw std::logic_error(std::string("Cannot calculate dot(mat[")
                             + std::to_string(sz[0])
                             + std::string(",")
                             + std::to_string(sz[1])
                             + std::string("], vec[")
                             + std::to_string(vec.size())
                             + std::string("])! Incorrect matrix or vector shapes!"));
    }
}

template <Vector V, Matrix M>
requires HaveGeneralType<V,M>
void check_dot_operands_sizes(const V& vec, const M& mat) {
    auto sz = size(mat);
    if (vec.size() != sz[0]) {
        throw std::logic_error(std::string("Cannot calculate dot(vec[")
                             + std::to_string(vec.size())
                             + std::string("], mat[")
                             + std::to_string(sz[0])
                             + std::string(",")
                             + std::to_string(sz[1])
                             + std::string("])! Incorrect vector or matrix shapes!"));
    }
}

template <Vector V1, ArithmeticVectorsLike<V1> V2>
void check_dot_operands_sizes(const V1& vec1, const V2& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::logic_error(std::string("Cannot calculate dot(vec[")
                             + std::to_string(vec1.size())
                             + std::string("], vec[")
                             + std::to_string(vec2.size())
                             + std::string("])! Incorrect vectors shapes!"));
    }
}

template <Matrix M1, ArithmeticVectorsLike<M1> M2>
void check_dotT_operands_sizes(const M1& mat1, const M2& mat2) {
    auto sz1 = size(mat1);
    auto sz2 = size(mat2);
    if (sz1[1] != sz2[1]) {
        throw std::logic_error(std::string("Cannot calculate dotT(mat[")
                             + std::to_string(sz1[0])
                             + std::string(",")
                             + std::to_string(sz1[1])
                             + std::string("], mat[")
                             + std::to_string(sz2[0])
                             + std::string(",")
                             + std::to_string(sz2[1])
                             + std::string("])! Incorrect matrices shapes!"));
    }
}

} // namespace detail

/* dot product of matrix and matrix TODO*/
template <Matrix M1, ArithmeticVectorsLike<M1> M2, ArithmeticVectorsLike<M1> M3>
void dot(const M1& mat1, const M2& mat2, M3& res) {
    detail::check_dot_operands_sizes(mat1,mat2);
    auto itres_row = res.begin();
    decltype(res[0].begin()) itres_col;

    auto itmat1_row     = mat1.begin()
        ,itmat1_row_end = mat1.end();
    decltype(mat1[0].begin()) itmat1_col;
    typename M1::basic_value_type m1_elem;
    
    decltype(mat2.begin()) itmat2_row, itmat2_row_end = mat2.end();
    decltype(mat2[0].begin()) itmat2_col, itmat2_col_end;
    
    while (itmat1_row < itmat1_row_end) {       // i
        itmat1_col = itmat1_row->begin();
        itmat2_row = mat2.begin();
        while (itmat2_row < itmat2_row_end) {   // k
            m1_elem = *itmat1_col;
            itres_col = itres_row->begin();
            itmat2_col = itmat2_row->begin();
            itmat2_col_end = itmat2_row->end();
            while(itmat2_col < itmat2_col_end) {                           // j
                *itres_col += m1_elem * (*itmat2_col);
                ++itmat2_col;
                ++itres_col;
            }
            ++itmat2_row;
            ++itmat1_col;
        }

        ++itmat1_row;
        ++itres_row;
    }
}

template <Matrix M1, ArithmeticVectorsLike<M1> M2>
auto dot(const M1& mat1, const M2& mat2)
        -> general_vector_type_t<M1,M2> {
    auto res = zeros<general_type_t<M1,M2>>(size(mat1)[0],size(mat2)[1]);
    dot(mat1,mat2,res);
    return res;
}


template <Matrix M, Vector V1, Vector V2>
requires HaveGeneralType<M,V1>
&& HaveGeneralType<V1,V2>
void dot(const M& mat, const V1& vec, V2& res) {
    detail::check_dot_operands_sizes(mat,vec);
    
    auto itres = res.begin();
    typename V2::basic_value_type temp_res;
    
    auto itmat_row = mat.begin()
        ,itmat_row_end = mat.end();
    decltype(mat[0].begin()) itmat_col;

    decltype(vec.begin()) itvec, itvec_end = vec.end();
    
    while (itmat_row < itmat_row_end) {
        temp_res = 0;
        itvec = vec.begin();
        itmat_col = itmat_row->begin();
        while (itvec < itvec_end) {
            temp_res += (*itmat_col) * (*itvec);
            ++itvec;
            ++itmat_col;
        }
        *itres = temp_res;
        ++itmat_row;
        ++itres;
    }

}


template <Matrix M, Vector V1, Vector V2>
requires HaveGeneralType<M,V1>
&& HaveGeneralType<V1,V2>
void dotInc(const M& mat, const V1& vec, V2& res) {
    detail::check_dot_operands_sizes(mat,vec);
    
    auto itres = res.begin();
    typename V2::basic_value_type temp_res;
    
    auto itmat_row = mat.begin()
        ,itmat_row_end = mat.end();
    decltype(mat[0].begin()) itmat_col;

    decltype(vec.begin()) itvec, itvec_end = vec.end();
    
    while (itmat_row < itmat_row_end) {
        temp_res = 0;
        itvec = vec.begin();
        itmat_col = itmat_row->begin();
        while (itvec < itvec_end) {
            temp_res += (*itmat_col) * (*itvec);
            ++itvec;
            ++itmat_col;
        }
        *itres += temp_res;
        ++itmat_row;
        ++itres;
    }

}


template <Matrix M, Vector V>
requires HaveGeneralType<M,V>
auto dot(const M& mat, const V& vec)
        -> vector_t<general_type_t<M,V>,1> {
    auto res = zeros<general_type_t<M,V>>(size(mat)[0]);
    dot(mat,vec,res);
    return res;
}

template <Vector V1, Matrix M, Vector V2>
requires HaveGeneralType<V1,M>
&& HaveGeneralType<V1,V2>
void dot(const V1& vec, const M& mat, V2& res) {
    detail::check_dot_operands_sizes(vec,mat);

    decltype(res.begin()) itres;

    auto itvec = vec.begin()
        ,itvec_end = vec.end();
    typename V1::basic_value_type temp_vec;

    auto itmat_row = mat.begin();
    decltype(mat[0].begin()) itmat_col, itmat_col_end;
    
    while (itvec < itvec_end) {
        itres = res.begin();
        itmat_col = itmat_row->begin();
        itmat_col_end = itmat_row->end();
        temp_vec = *itvec;
        while (itmat_col < itmat_col_end) {
            *itres += (*itmat_col)*temp_vec;
            ++itmat_col;
            ++itres;
        }

        ++itvec;
        ++itmat_row;
    }
}

template <Matrix M, Vector V>
requires HaveGeneralType<M,V>
auto dot(const V& vec, const M& mat)
        -> vector_t<general_type_t<M,V>,1> {
    auto res = zeros<general_type_t<M,V>>(size(mat)[1]);
    dot(vec,mat,res);
    return res;
}


template <Vector V1, ArithmeticVectorsLike<V1> V2>
auto dot(const V1& vec1, const V2& vec2)
        -> general_type_t<V1,V2> {
    detail::check_dot_operands_sizes(vec1,vec2);
    general_type_t<V1,V2> res{};

    auto itvec1 = vec1.begin()
        ,itvec1_end = vec1.end();
    auto itvec2 = vec2.begin();

    while (itvec1 < itvec1_end) {
        res += (*itvec1) * (*itvec2);
        ++itvec1;
        ++itvec2;
    }
    return res;
}

template <Matrix M1, ArithmeticVectorsLike<M1> M2, ArithmeticVectorsLike<M1> M3>
void dotT(const M1& mat1, const M2& mat2T, M3& res) {
    detail::check_dotT_operands_sizes(mat1,mat2T);

    auto itres_row = res.begin();
    decltype(res[0].begin()) itres_col;
    typename M3::basic_value_type temp_res;

    auto itmat1_row = mat1.begin()
        ,itmat1_row_end = mat1.end();
    decltype(mat1[0].begin()) itmat1_col, itmat1_col_end;

    decltype(mat2T.begin()) itmat2_row
                           ,itmat2_row_end = mat2T.end();
    decltype(mat2T[0].begin()) itmat2_col;
    
    while (itmat1_row < itmat1_row_end) {
        itmat2_row = mat2T.begin();
        itres_col = itres_row->begin();
        while (itmat2_row < itmat2_row_end) {
            temp_res = 0;
            itmat1_col = itmat1_row->begin();
            itmat1_col_end = itmat1_row->end();
            itmat2_col = itmat2_row->begin();
            while (itmat1_col < itmat1_col_end) {
                temp_res += (*itmat1_col) * (*itmat2_col);
                ++itmat1_col;
                ++itmat2_col;
            }
            *itres_col = temp_res;
            ++itmat2_row;
            ++itres_col;
        }

        ++itmat1_row;
        ++itres_row;
    }
}

template <Matrix M1, ArithmeticVectorsLike<M1> M2>
auto dotT(const M1& mat1, const M2& mat2T)
        -> general_vector_type_t<M1,M2> {
    auto res = zeros<general_type_t<M1,M2>>(size(mat1)[0],size(mat2T)[0]);
    dotT(mat1,mat2T,res);
    return res;
}

template <Matrix M1, ArithmeticVectorsLike<M1> M2, ArithmeticVectorsLike<M1> M3>
void dotUL(const M1& U, const M2& L, M3& res) {
    detail::check_dotT_operands_sizes(U,L);
    typename M3::basic_value_type tmp;
    auto n = U.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            tmp = 0;
            for (size_t k = i<j?j:i; k < n; ++k) {
                tmp += U[i][k] * L[k][j];
            }
            res[i][j] = tmp;
        }
    }
}

template <Matrix M1, ArithmeticVectorsLike<M1> M2>
auto dotUL(const M1& U, const M2& L)
        -> general_vector_type_t<M1,M2> {
    auto res = zeros<general_type_t<M1,M2>>(size(U)[0],size(L)[0]);
    dotUL(U,L,res);
    return res;
}

template <Vector V, Matrix M, VectorLike V3D>
requires HaveGeneralType<V,M>
&& HaveGeneralType<V,V3D>
&& Matrix<typename V3D::value_type>
void dot(const V3D& v3d
       , const V& v
             , M& res)
{
    detail::check_dot_operands_sizes(v3d,v);

    decltype(res.begin()) itres_row = res.begin();
    decltype(res[0].begin()) itres_col;
    typename M::basic_value_type res_temp;

    decltype(v.begin()) itv, itv_end = v.end();

    auto itM0 = v3d.begin(), itM0_end = v3d.end();
    decltype(v3d[0].begin()) itM1, itM1_end;
    decltype(v3d[0][0].begin()) itM2;


    while (itM0 < itM0_end) {
        itM1 = itM0->begin();
        itM1_end = itM0->end();
        itres_col = itres_row->begin();
        while (itM1 < itM1_end) {
            res_temp = 0.;
            itM2 = itM1->begin();
            itv = v.begin();
            while (itv < itv_end) {
                res_temp += (*itM2) * (*itv);
                ++itM2;
                ++itv;
            }
            *itres_col = res_temp;
            ++itM1;
            ++itres_col;
        }
    
        ++itM0;
        ++itres_row;
    }
}

template <Vector V, Matrix M, VectorLike V3D>
requires HaveGeneralType<V,M>
&& HaveGeneralType<V,V3D>
&& Matrix<typename V3D::value_type>
void dot(const V& v, const V3D& v3d, M& res)
{
    detail::check_dot_operands_sizes(v,v3d);

    decltype(res.begin()) itres_row;
    decltype(res[0].begin()) itres_col;

    auto itv = v.begin();
    typename V::basic_value_type vk;

    auto itM0 = v3d.begin(), itM0_end = v3d.end();
    decltype(v3d[0].begin()) itM1, itM1_end;
    decltype(v3d[0][0].begin()) itM2, itM2_end;


    while (itM0 < itM0_end) {
        vk = *itv;
        itM1 = itM0->begin();
        itM1_end = itM0->end();
        itres_row = res.begin();
        while (itM1 < itM1_end) {
            itM2 = itM1->begin();
            itM2_end = itM1->end();
            itres_col = itres_row->begin();
            while (itM2 < itM2_end) {
                *itres_col += vk * (*itM2);
                ++itM2;
                ++itres_col;
            }
            ++itM1;
            ++itres_row;
        }
    
        ++itM0;
        ++itv;
    }
    /*  for (size_t k = 0) {
            for (size_t i = 0) {
                for (size_t j = 0) {
                    res[i][j] += v3d[k][i][j] * v[k]
                }
            }
        }
    */
}

template <Vector V, VectorLike V3D>
requires HaveGeneralType<V,V3D>
&& Matrix<typename V3D::value_type>
auto dot(const V3D& v3d, const V& v)
        -> vector_t<general_type_t<V,V3D>,2> 
{
    vector_t<general_type_t<V,V3D>,2> res = zeros<general_type_t<V,V3D>>(v3d.size(),v3d[0].size());
    dot(v3d,v,res);
    return res;
}

template <Vector V, VectorLike V3D>
requires HaveGeneralType<V,V3D>
&& Matrix<typename V3D::value_type>
auto dot(const V& v, const V3D& v3d)
        -> vector_t<general_type_t<V,V3D>,2> 
{
    vector_t<general_type_t<V,V3D>,2> res = zeros<general_type_t<V,V3D>>(v3d[0].size(),v3d[0][0].size());
    dot(v,v3d,res);
    return res;
}


template <Vector V1, ArithmeticVectorsLike<V1> V2>
auto cross(const V1& vec1, const V2& vec2)
        -> general_vector_type_t<V1,V2> {
    general_vector_type_t<V1,V2> res = {
        vec1[1]*vec2[2] - vec1[2]*vec2[1],
        vec1[2]*vec2[0] - vec1[0]*vec2[2],
        vec1[0]*vec2[1] - vec1[1]*vec2[0]
    };
    return res;
}



template <Vector V1, ArithmeticVectorsLike<V1> V2>
auto kron(const V1& v1, const V2& v2)
        -> vector_t<general_type_t<V1,V2>,1>
{
    vector_t<general_type_t<V1,V2>,1> res = zeros<general_type_t<V1,V2>>(v1.size()*v2.size());
    typename V1::basic_value_type v1_value;

    auto it =       res.begin();
   
    auto it1 =      v1.begin(),     end1 = v1.end();
   
    auto end2 =     v2.end();
    decltype(end2)  it2;

    while (it1 != end1) {
        v1_value = *it1;
        it2 = v2.begin();
        while (it2 != end2) {
            *it = v1_value * (*it2);
            ++it;   ++it2;
        }
        ++it1;
    }
    return res;
}

template <Vector V, Matrix M>
requires HaveGeneralType<V,M>
auto kron(const V& v, const M& m)
        -> vector_t<general_type_t<V,M>,2>
{
    auto m_sz = size(m);
    vector_t<general_type_t<V,M>,2> res = zeros<general_type_t<V,M>>(v.size()*m_sz[0], m_sz[1]);
    typename V::basic_value_type v_value;

    auto it_row = res.begin();
    decltype(res[0].begin()) it_col; 
    
    auto itv = v.begin(), endv = v.end();
    
    auto endm_row = m.end();
    decltype(endm_row) itm_row;
    decltype(m[0].end()) itm_col, endm_col;

    while (itv != endv) {

        v_value = *itv;
        itm_row = m.begin();

        while (itm_row != endm_row) {

            itm_col = (*itm_row).begin();
            endm_col = (*itm_row).end();
            it_col = (*it_row).begin();

            while (itm_col != endm_col) {

                *it_col = v_value * (*itm_col);
                ++it_col;   ++itm_col;
            
            }
            ++it_row;   ++itm_row;

        }
        ++itv;
    }
    return res;
}

template <Matrix M1, ArithmeticVectorsLike<M1> M2, ArithmeticVectorsLike<M1> M3>
requires (!std::is_const_v<M3>)
void kron(const M1& m1, const M2& m2, M3& res) {
    auto m1_sz = size(m1);
    auto m2_sz = size(m2);
    typename M1::basic_value_type m1_value;

    decltype(res.begin())    it_row;
    decltype(res[0].begin()) it_col;

    decltype(m2.begin())    it2_row, end2_row = m2.end();
    decltype(m2[0].begin()) it2_col, end2_col;  
    
    size_t row_m1 = 0, col_m1 = 0;
    
    for (row_m1 = 0; row_m1 < m1_sz[0]; ++row_m1) {

        for (col_m1 = 0; col_m1 < m1_sz[1]; ++col_m1) {
            
            m1_value = m1[row_m1][col_m1];
            it_row = res.begin()+row_m1*m2_sz[0];
            it2_row = m2.begin();
            
            while (it2_row != end2_row) {

                it_col = (*it_row).begin() + col_m1*m2_sz[1];
                it2_col = (*it2_row).begin();
                end2_col = (*it2_row).end();
                
                while (it2_col != end2_col) {

                    *it_col = m1_value * (*it2_col);
                    ++it_col;   ++it2_col;
                
                }
                ++it_row;   ++it2_row;
            }
            
        }
        
    }

}

template <Matrix M1, ArithmeticVectorsLike<M1> M2>
auto kron(const M1& m1, const M2& m2)
        -> vector_t<general_type_t<M1,M2>,2>
{
    auto m1_sz = size(m1);
    auto m2_sz = size(m2);
    vector_t<general_type_t<M1,M2>,2> res = 
                                zeros<general_type_t<M1,M2>>(m1_sz[0]*m2_sz[0],
                                                             m1_sz[1]*m2_sz[1]);
    kron(m1,m2,res);
    return res;
}


/* addkron(m1,m2,res) is the same as res += kron(m1,m2) */
template <Matrix M1, ArithmeticVectorsLike<M1> M2, ArithmeticVectorsLike<M1> M3>
requires (!std::is_const_v<M3>)
void addkron(const M1& m1, const M2& m2, M3& res) {
    auto m1_sz = size(m1);
    auto m2_sz = size(m2);
    typename M1::basic_value_type m1_value;

    decltype(res.begin())    it_row;
    decltype(res[0].begin()) it_col;

    decltype(m2.begin())    it2_row, end2_row = m2.end();
    decltype(m2[0].begin()) it2_col, end2_col;  
    
    size_t row_m1 = 0, col_m1 = 0;
    
    for (row_m1 = 0; row_m1 < m1_sz[0]; ++row_m1) {

        for (col_m1 = 0; col_m1 < m1_sz[1]; ++col_m1) {
            
            m1_value = m1[row_m1][col_m1];
            it_row = res.begin()+row_m1*m2_sz[0];
            it2_row = m2.begin();
            
            while (it2_row != end2_row) {

                it_col = (*it_row).begin() + col_m1*m2_sz[1];
                it2_col = (*it2_row).begin();
                end2_col = (*it2_row).end();
                
                while (it2_col != end2_col) {

                    *it_col += m1_value * (*it2_col);
                    ++it_col;   ++it2_col;
                
                }
                ++it_row;   ++it2_row;
            }
            
        }
        
    }

}


template <NotVectorLike T>
auto cos(const T& value)
        -> decltype(std::cos(value))
{
    return std::cos(value);
}

template <VectorLike V>
auto cos(const V& vec)
        -> std::decay_t<V>
{
    std::decay_t<V> res = vec;
    auto it = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it < end) {
        *it = cos(*it1);
        ++it; ++it1;
    }
    return res;
}

template <VectorLike V>
requires (!std::is_lvalue_reference_v<V>)
auto cos(V&& vec)
        -> std::remove_reference_t<V>
{
    std::remove_reference_t<V> res = std::move(vec);
    auto it = res.begin(), end = res.end();
    while (it < end) {
        *it = cos(*it);
        ++it;
    }
    return res;
}



template <NotVectorLike T>
auto sin(const T& value)
        -> decltype(std::sin(value))
{
    return std::sin(value);
}

template <VectorLike V>
auto sin(const V& vec)
        -> std::decay_t<V>
{
    std::decay_t<V> res = vec;
    auto it = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it < end) {
        *it = sin(*it1);
        ++it; ++it1;
    }
    return res;
}

template <VectorLike V>
requires (!std::is_lvalue_reference_v<V>)
auto sin(V&& vec)
        -> std::remove_reference_t<V>
{
    std::remove_reference_t<V> res = std::move(vec);
    auto it = res.begin(), end = res.end();
    while (it < end) {
        *it = sin(*it);
        ++it;
    }
    return res;
}

template <NotVectorLike T>
auto tan(const T& value)
        -> decltype(std::tan(value))
{
    return std::tan(value);
}

template <VectorLike V>
auto tan(const V& vec)
        -> std::decay_t<V>
{
    std::decay_t<V> res = vec;
    auto it = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it < end) {
        *it = tan(*it1);
        ++it; ++it1;
    }
    return res;
}

template <VectorLike V>
requires (!std::is_lvalue_reference_v<V>)
auto tan(V&& vec)
        -> std::remove_reference_t<V>
{
    std::remove_reference_t<V> res = std::move(vec);
    auto it = res.begin(), end = res.end();
    while (it < end) {
        *it = tan(*it);
        ++it;
    }
    return res;
}



template <NotVectorLike T>
vector_t<T, 1> linspace(T from, T to, size_t num) {
    auto res = zeros<T>(num);
    T cur = 0;
    T num_ = num;
    auto it = res.begin();
    for (; cur < num_; ++cur, ++it) {
        *it = from + (to - from)*cur/(num_-1);
    }
    return res;
}



template <std::input_iterator InputIt>
requires NotVectorLike<typename InputIt::value_type>
auto sumsq(InputIt first, InputIt last)
		-> typename InputIt::value_type {
    typename InputIt::value_type res = 0;
    while (first < last) {
        res += (*first) * (*first);
        ++first;
    }
    return res;
}
template <std::input_iterator InputIt>
requires VectorLike<typename InputIt::value_type>
auto sumsq(InputIt first, InputIt last)
		-> basic_value_type_t<typename InputIt::value_type>
{
    basic_value_type_t<typename InputIt::value_type> res = 0;
    while (first < last) {
        res += sumsq(first->begin(), first->end());
        ++first;
    }
    return res;
}

template <std::input_iterator InputIt>
requires NotVectorLike<typename InputIt::value_type>
double norm(InputIt first, InputIt last) {
	return std::sqrt(sumsq(first,last));
}


template <std::input_iterator InputIt>
requires VectorLike<typename InputIt::value_type>
double norm(InputIt first, InputIt last) {
	return std::sqrt(sumsq(first,last));
}


template <VectorLike V>
typename V::basic_value_type sumsq(const V& vec) {
    return sumsq(vec.begin(),vec.end());
}

template <VectorLike V>
double norm(const V& vec) {
    return norm(vec.begin(),vec.end());
}





template <NotVectorLike T>
auto pow(const T& value, double n)
        -> decltype(std::pow(value,n))
{
    return std::pow(value,n);
}

template <VectorLike V>
auto pow(const V& vec, double n)
        -> std::decay_t<V>
{
    std::decay_t<V> res = vec;
    auto it = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it < end) {
        *it = pow(*it1,n);
        ++it; ++it1;
    }
    return res;
}

template <VectorLike V>
requires (!std::is_lvalue_reference_v<V>)
auto pow(V&& vec, double n)
        -> std::remove_reference_t<V>
{
    std::remove_reference_t<V> res = std::move(vec);
    auto it = res.begin(), end = res.end();
    while (it < end) {
        *it = pow(*it,n);
        ++it;
    }
    return res;
}


template <NotVectorLike T>
auto abs(const T& value)
        -> decltype(std::abs(value))
{
    return std::abs(value);
}

template <VectorLike V>
auto abs(const V& vec)
        -> std::decay_t<V>
{
    std::decay_t<V> res = vec;
    auto it = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it < end) {
        *it = abs(*it1);
        ++it; ++it1;
    }
    return res;
}

template <VectorLike V>
requires (!std::is_lvalue_reference_v<V>)
auto abs(V&& vec)
        -> std::remove_reference_t<V>
{
    std::remove_reference_t<V> res = std::move(vec);
    auto it = res.begin(), end = res.end();
    while (it < end) {
        *it = abs(*it);
        ++it;
    }
    return res;
}

/*Linear solver using Cholesky decomposition: for symmetrix and positive defined matrices*/
void solve_llt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x, vector_t<double,2>& L, size_t n);
void solve_llt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x, size_t n);
void solve_llt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x, vector_t<double,2>& L);
void solve_llt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x);
// // for band-matrices
// void solve_llt(const vector_t<double,1>& A, const vector<size_t>& diag, vector<double>& b, vector<double>& x, size_t n, vector_t<double,2>& L);
// void solve_llt(const vector_t<double,1>& A, const vector<size_t>& diag, vector<double>& b, vector<double>& x, size_t n);
// void solve_llt(const vector_t<double,1>& A, const vector<size_t>& diag, vector<double>& b, vector<double>& x, vector_t<double,2>& L);
// void solve_llt(const vector_t<double,1>& A, const vector<size_t>& diag, vector<double>& b, vector<double>& x);

/*Cholesky decomposition: for symmetrix and positive defined matrices*/
void LLT(const vector_t<double,2>& A, vector_t<double,2>& L, size_t n);
void LLT(const vector_t<double,2>& A, vector_t<double,2>& L);
// // for band-matrices
// void LLT(const vector_t<double,1>& A, const vector<size_t>& diag, vector<double>& L, size_t n);
// void LLT(const vector_t<double,1>& A, const vector<size_t>& diag, vector<double>& L);



/*Linear solver using LDLT decomposition: for symmetrix matrices*/
void solve_ldlt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x, vector_t<double,2>& L, vector<double>& D, size_t n);
void solve_ldlt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x);

void solve_ldlt(const vector_t<double,1>& A, const vector<size_t>& diag, const vector<double>& b, vector<double>& x, vector<double>& LT, vector<double>& D, size_t n);
void solve_ldlt(const vector_t<double,1>& A, const vector<size_t>& diag, const vector<double>& b, vector<double>& x);


void LDLT(const vector_t<double,2>& A, vector_t<double,2>& L, vector<double>& D, vector<double>& g, size_t n);
void LDLT(const vector_t<double,2>& A, vector_t<double,2>& L, vector<double>& D);

void LDLT(const vector_t<double,1>& A, const vector<size_t>& diag, vector<double>& LT, vector<double>& D, vector<double>& g, size_t n);
void LDLT(const vector_t<double,1>& A, const vector<size_t>& diag, vector<double>& LT, vector<double>& D);


/*Linear solver using LU decomposition*/
void solve_lu(vector_t<double,2>& A, const vector<double>& b, vector<double>& x, size_t n);
void solve_lu(vector_t<double,2>& A, const vector<double>& b, vector<double>& x);

void LU(vector_t<double,2>& A, size_t n);
void LU(vector_t<double,2>& A);


/*Linear solver using pseudo-inverse matrix. Find solution with minimal norm*/
void psolve(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x, vector_t<double,2>& AAT, vector_t<double,2>& L, vector<double>& D);
void psolve(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x);

// void psolve_weight(const vector_t<double,2>& A, const vector<double>& b, const vector_t<double,2>& W, vector<double>& x, vector_t<double,2>& AAT, vector_t<double,2>& L, vector<double>& D);
void psolve_weight(const vector_t<double,2>& A, const vector<double>& b, const vector_t<double,2>& invW, vector<double>& x);


/*Transpose rectangular matrix*/
template <Matrix M1, Matrix M2>
requires (!std::is_const_v<M2>)
void transpose(const M1& A, M2& AT) {
    auto sz_A = size(A);
    if (sz_A[0] == sz_A[1]) {
        transpose_square(A,AT);
        return;
    }

    for (size_t i = 0; i < sz_A[0]; ++i) {
        for (size_t j = 0; j < sz_A[1]; ++j) {
            AT[j][i] = A[i][j];
        }
    }
}

template <Matrix M>
auto transpose(const M& A)
        -> vector_t<typename M::basic_value_type,2>
{
    vector_t<typename M::basic_value_type,2> AT = zeros<typename M::basic_value_type>(A[0].size(),A.size());
    transpose(A,AT);
    return AT;
}

template <Matrix M1, Matrix M2>
requires (!std::is_const_v<M2>)
void transpose_square(const M1& A, M2& AT) {
    for (size_t i = 0; i < A.size(); ++i) {
        AT[i][i] = A[i][i];
        for (size_t j = i+1; j < A.size(); ++j) {
            AT[i][j] = A[j][i];
            AT[j][i] = A[i][j];
        }
    }
}

template <Matrix M>
auto transpose_square(const M& A)
        -> vector_t<typename M::basic_value_type,2>
{
    vector_t<typename M::basic_value_type,2> AT = zeros<typename M::basic_value_type>(A.size(),A.size());
    transpose_square(A,AT);
    return AT;
}

void invLowTri(const vector_t<double,2>& A, vector_t<double,2>& invA);
void invLowTriUnit(const vector_t<double,2>& A, vector_t<double,2>& invA);
void invUpTri(const vector_t<double,2>& A, vector_t<double,2>& invA);
void inv(const vector_t<double,2>& A, vector_t<double,2>& invA);
void invSym(const vector_t<double,2>& A, vector_t<double,2>& invA
            , vector_t<double,2>& L, vector<double>& D
            , vector<double>& temp, vector_t<double,2>& invL, size_t n);
void invSym(const vector_t<double,2>& A, vector_t<double,2>& invA);
/*Inverse rectangular matrix: pseudo-inverse matrix*/
void pinv(const vector_t<double,2>& A, vector_t<double,2>& pinvA);




template <Matrix M>
auto matrix_block(const M& mat, size_t i, size_t j, size_t block_rows, size_t block_cols)
		-> math::vector< Slice<decltype(mat[0].begin()),
							   decltype(mat[0].begin())> > // const_block_t
{
	math::vector< Slice<decltype(mat[0].begin()),
						decltype(mat[0].begin())> > block(block_rows);
	for (size_t k = 0; k < block_rows; ++k) {
		block[k].new_slice(mat[block_rows*i+k].begin()+block_cols*j,
				           mat[block_rows*i+k].begin()+block_cols*(j+1) );
	}
	return block;
}

template <Matrix M>
auto matrix_block(M& mat, size_t i, size_t j, size_t block_rows, size_t block_cols)
		-> math::vector< Slice<decltype(mat[0].begin()),
							   decltype(mat[0].begin())> > // block_t
{
	math::vector< Slice<decltype(mat[0].begin()),
						decltype(mat[0].begin())> > block(block_rows);
	for (size_t k = 0; k < block_rows; ++k) {
		block[k].new_slice(mat[block_rows*i+k].begin()+block_cols*j,
				           mat[block_rows*i+k].begin()+block_cols*(j+1) );
	}
	return block;
}

template <Matrix M>
void matrix_block_set(M& mat, size_t i, size_t j, size_t block_rows, size_t block_cols
					,math::vector< Slice<decltype(mat[0].begin()),
									     decltype(mat[0].begin())> >& block) // block_t
{
	for (size_t k = 0; k < block_rows; ++k) {
		block[k].new_slice(mat[block_rows*i+k].begin()+block_cols*j,
				           mat[block_rows*i+k].begin()+block_cols*(j+1) );
	}
}


template <Vector V1, ArithmeticVectorsLike<V1> V2>
void assignment_sub(const V1& source, V2& dest) {
    auto s = source.begin()
        ,send =source.end();
    auto d = dest.begin();

    while(s < send) {
        *d = -(*s);
        ++s;
        ++d;
    }
}

template <VectorLike V1, ArithmeticVectorsLike<V1> V2>
void assignment_sub(const V1& source, V2& dest) {
    auto s = source.begin()
        ,send =source.end();
    auto d = dest.begin();

    while(s < send) {
        assignment_sub(*s,*d);
        ++s;
        ++d;
    }
}

       
} // namespace math 