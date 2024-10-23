#pragma once

#include "vector.hpp"
#include "slice.hpp"
#include <cmath>
#include <numbers>
#include <string>

namespace math {

constexpr double pi = std::numbers::pi;

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

double det_(const vector_t<double,2>& mat);

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

template <Matrix M1, ArithmeticVectorsLike<M1> M2>
auto dot(const M1& mat1, const M2& mat2)
        -> general_vector_type_t<M1,M2> {
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
        throw std::logic_error(std::string("Cannot calculate dot(mat[")
                             + std::to_string(sz_mat[0])
                             + std::string(",")
                             + std::to_string(sz_mat[1])
                             + std::string("], vec[")
                             + std::to_string(vec.size())
                             + std::string("])! Incorrect matrix or vector shape!"));
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

template <Matrix M, Vector V, Vector V0>
requires HaveGeneralType<M,V>
&& HaveGeneralType<V,V0>
void dot(const M& mat, const V& vec, V0& res) {
    auto sz_mat = size(mat);
    if (sz_mat[1] != vec.size()) {
        throw std::logic_error(std::string("Cannot calculate dot(mat[")
                             + std::to_string(sz_mat[0])
                             + std::string(",")
                             + std::to_string(sz_mat[1])
                             + std::string("], vec[")
                             + std::to_string(vec.size())
                             + std::string("])! Incorrect matrix or vector shape!"));
    }
    auto it_res = res.begin();
    decltype(vec.begin()) it_vec;
    auto it_mat_row = mat.begin(), it_mat_row_end = mat.end();
    decltype(mat[0].begin()) it_mat_col,it_mat_col_end;

    typename V0::basic_value_type tmp;
    while (it_mat_row != it_mat_row_end) {
        it_vec = vec.begin();
        it_mat_col = it_mat_row->begin();
        it_mat_col_end = it_mat_row->end();
        tmp = 0;
        while (it_mat_col != it_mat_col_end) {
            tmp += (*it_mat_col) * (*it_vec);
            ++it_mat_col; ++it_vec;
        }
        *it_res = tmp;
        ++it_mat_row; ++it_res;
    }
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

template <Matrix M1, ArithmeticVectorsLike<M1> M2>
auto kron(const M1& m1, const M2& m2)
        -> vector_t<general_type_t<M1,M2>,2>
{
    auto m1_sz = size(m1);
    auto m2_sz = size(m2);
    vector_t<general_type_t<M1,M2>,2> res = zeros<general_type_t<M1,M2>>(m1_sz[0]*m2_sz[0],
                                                                         m1_sz[1]*m2_sz[1]);
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

    return res;
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
    size_t cur = 0;
    auto it = res.begin();
    for (; cur < num; ++cur, ++it) {
        *it = from + (to - from)*cur/(num-1);
    }
    return res;
}



template <std::input_iterator InputIt>
requires NotVectorLike<typename InputIt::value_type>
auto sumsq(InputIt first, InputIt last)
		-> typename InputIt::value_type {
    typename InputIt::value_type res = 0;
    while (first < last) {
        res += *first * *first;
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

vector<double> solve(const vector_t<double,2>& A, const vector<double>& b);
vector<double> solve2(const vector_t<double,2>& A, const vector<double>& b);
vector<double> solve2(const vector<double>& A, const vector<double>& b, const vector<size_t>& diag);
vector<double> psolve(const vector_t<double,2>& A, const vector<double>& b);

void LDLT(const vector_t<double,2>& K, vector<double>& D, vector_t<double,2>& L);
void LU(const vector_t<double,2>& K, vector_t<double,2>& L, vector_t<double,2>& U);
void gauss_backward_triup(vector<double>& x, const vector_t<double,2>& A, const vector<double>& b);
void gauss_backward_tridown(vector<double>& x, const vector_t<double,2>& A, const vector<double>& b);

vector_t<double,2> inv(const vector_t<double,2>& A);
vector_t<double,2> pinv(const vector_t<double,2>& A);

template <Matrix M>
auto transpose(const M& A)
        -> vector_t<typename M::basic_value_type,2>
{
    auto sz_A = size(A);
    if (sz_A[0] == sz_A[1])
        return transpose_square(A);

    vector_t<typename M::basic_value_type,2> AT = zeros<typename M::basic_value_type>(sz_A[1],sz_A[0]);
    for (size_t i = 0; i < sz_A[0]; ++i) {
        for (size_t j = 0; j < sz_A[1]; ++j) {
            AT[j][i] = A[i][j];
        }
    }
    return AT;
}
template <Matrix M>
auto transpose_square(const M& A)
        -> vector_t<typename M::basic_value_type,2>
{
    vector_t<typename M::basic_value_type,2> AT = zeros<typename M::basic_value_type>(A.size(),A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        AT[i][i] = A[i][i];
        for (size_t j = i+1; j < A.size(); ++j) {
            AT[i][j] = A[j][i];
            AT[j][i] = A[i][j];
        }
    }
    return AT;
}

vector_t<double,2> invLowTri(const vector_t<double,2>& A);
vector_t<double,2> invUpTri(const vector_t<double,2>& A);
vector_t<double,2> inv(const vector_t<double,2>& A);


vector_t<double,2> rotation_tensor(double theta);

/* Vector invariant of the dot product of 2 tensors (square matrises, size=3):
 mat1 and transpose(mat2) */
template <Matrix M1, ArithmeticVectorsLike<M1> M2>
auto vector_invariant(const M1& mat1, const M2& mat2_toT)
        -> vector<general_type_t<M1,M2>>
{
    
    vector<general_type_t<M1,M2>> res = {
		math::dot(mat1[1],mat2_toT[2]) - math::dot(mat1[2],mat2_toT[1]),
        math::dot(mat1[2],mat2_toT[0]) - math::dot(mat1[0],mat2_toT[2]),
        math::dot(mat1[0],mat2_toT[1]) - math::dot(mat1[1],mat2_toT[0]),
	};
	return res;
}


template <Matrix M>
auto matrix_block(const M& mat, size_t i, size_t j, size_t block_rows, size_t block_cols)
		-> math::vector< Slice<decltype(mat[0].begin()),
							   decltype(mat[0].begin())> >
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
							   decltype(mat[0].begin())> >
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
									     decltype(mat[0].begin())> >& block)
{
	for (size_t k = 0; k < block_rows; ++k) {
		block[k].new_slice(mat[block_rows*i+k].begin()+block_cols*j,
				           mat[block_rows*i+k].begin()+block_cols*(j+1) );
	}
}

template <Vector V>
auto skew_symmetric_tensor(const V& vec)
		-> vector_t<typename V::basic_value_type,2>
{
	vector_t<typename V::basic_value_type,2> res = {
		{   0   ,-vec[2], vec[1]},
		{ vec[2],   0   ,-vec[0]},
		{-vec[1], vec[0],   0   }
	};
	return res;
}
} // namespace math 