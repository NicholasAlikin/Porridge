#include "vector.h"

namespace math {

/*==================================
// Generate vector-like objects
==================================*/

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

template <typename T, typename... Sizes>
requires (std::is_convertible_v<Sizes,size_t> && ...)
auto zeros(Sizes... sizes) 
        -> vector_t<T,sizeof...(Sizes)> {
    return zeros_helper<T,sizeof...(Sizes)>(sizes...);
}

template <typename T, size_t... sizes>
auto zeros(std::integer_sequence<size_t,sizes...> int_seq) 
        -> vector_t<T,sizeof...(sizes)> {
    return zeros_helper<T,int_seq.size()>((sizes, int_seq.size()) ...);
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

template <typename T, size_t Dim>
vector_t<T,Dim> zeros(const std::array<size_t,Dim>& sizes) {
    return zeros_helper<T,Dim>(sizes.begin());
}

template <typename T, VectorLike Vec>
auto zeros(const Vec& v)
        -> vector_t<T,vector_dim_v<Vec>> {
    return zeros<T>(size(v));
}

// repmat - vector-like full of given value with type vector::basic_value_type

template <typename T, size_t Dim>
requires (Dim == 1)
vector_t<T,Dim> repmat_helper(const T& value, size_t sz) {
    return vector_t<T,Dim>(sz,value);
}
template <typename T, size_t Dim>
vector_t<T,Dim> repmat_helper(const T& value, size_t sz, auto... sizes) {
    return vector_t<T,Dim>(sz, repmat_helper<T,Dim-1>(value, sizes...));
}

template <typename T, typename... Sizes>
requires (std::is_convertible_v<Sizes,size_t> && ...)
vector_t<T,sizeof...(Sizes)> repmat(const T& value, Sizes... sizes) {
    return repmat_helper<T,sizeof...(Sizes)>(value, sizes...);
}


template <typename T, size_t Dim>
requires (Dim == 1)
vector_t<T,Dim> repmat_helper(const T& value, const size_t* sizes) {
    return vector_t<T,Dim>(*sizes, value);
}

template <typename T, size_t Dim>
vector_t<T,Dim> repmat_helper(const T& value, const size_t* sizes) {
    const size_t sz = *sizes;
    return vector_t<T,Dim>(sz, repmat_helper<T,Dim-1>(value, ++sizes));
}

template <typename T, size_t Dim>
vector_t<T,Dim> repmat(const T& value, std::array<size_t,Dim> sizes) {
    return repmat_helper<T,Dim>(value, sizes.begin());
}

template <typename T, VectorLike Vec>
auto repmat(const T& value, const Vec& v)
        -> vector_t<T,vector_dim_v<Vec>> {
    return repmat<T>(value, size(v));
}

// eye - vector_t<T,2> equal to zeros except of diagonal elements
// which equal to 1
template <typename T, size_t sz>
vector_t<T,2> eye() {
    vector_t<T,2> v = zeros<T>(sz,sz);
    for (size_t i = 0; i < sz; ++i) {
        v[i][i] = 1;
    }
    return v;
}

template <typename T, size_t sz>
requires (sz == 1)
vector_t<T,1> eye() {
    return vector_t<T,1>(1,1);
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

template <Matrix M>
auto det(const M& mat)
        -> typename M::basic_value_type {
    auto sz = size(mat);
    if (sz[0] != sz[1]) {
        throw std::logic_error("Cannot calculate matrix determinant! Matrix dimensions must be equal!");
    }

    return det_(mat);
}


template <Matrix M1, ArithmeticVectorsLike<M1> M2>
auto dot(const M1& mat1, const M2& mat2)
        -> general_arithmetic_vector_like_type_t<M1,M2> {
    auto sz1 = size(mat1);
    auto sz2 = size(mat2);
    if (sz1[1] != sz2[0]) {
        throw std::logic_error("Cannot calculate dot(mat1,mat2)! Incorrect matrices shapes!");
    }
    auto res = zeros<general_type_t<M1,M2>>(sz1[0],sz2[1]);
    size_t i,j,k;
    for (i = 0; i < sz1[0]; ++i) {
        for (k = 0; k < sz1[1]; ++k) {
            for (j = 0; j < sz2[1]; ++j) {
                res[i][j] += mat1[i][k]*mat2[k][j];
            }
        }  
    }
    return res;
}

template <Matrix M, Vector V>
requires HaveGeneralType<M,V>
auto dot(const M& mat, const V& vec)
        -> vector_t<general_type_t<M,V>,1> {
    auto sz_mat = size(mat);
    if (sz_mat[1] != vec.size()) {
        throw std::logic_error("Cannot calculate dot(mat,vec)! Incorrect matrix or vector shape!");
    }
    auto res = zeros<general_type_t<M,V>>(sz_mat[0]);
    size_t i,j;
    for (i = 0; i < sz_mat[0]; ++i) {
        for (j = 0; j < sz_mat[1]; ++j) {
            res[i] += mat[i][j]*vec[j];
        }
    }
    return res;
}

template <Vector V1, ArithmeticVectorsLike<V1> V2>
auto dot(const V1& vec1, const V2& vec2)
        -> general_type_t<V1,V2> {
    if (vec1.size() != vec2.size()) {
        throw std::logic_error("Cannot calculate dot(vec1,vec2)! Incorrect vectors shapes!");
    }
    general_type_t<V1,V2> res{};
    for (size_t i = 0; i < vec1.size(); ++i) {
        res += vec1[i]*vec2[i];
    }
    return res;
}






} // namespace math 