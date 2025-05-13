#pragma once
/*Function specified for vector like objects:
    o operatos (+,-,/,*,>>)
    o copy (without allocations)
    o size
    */

#include "vector.hpp"

namespace math {

/*==================================
// vector mathematical operators overloadings
==================================*/
/*

Operators are performed componentwise.

There are two main groups of overloadings
based on the operands types:
    1 - both operands types are vector-like
    2 - one operand type is vector-like,
        second - basic_value_type.

Each group containts move-semantic supported overloadings.

*/

/*==================================
// operator+, +=
==================================*/

template <VectorLike Vec1, VectorsSameDim<Vec1> Vec2>
requires std::convertible_to<basic_value_type_t<Vec2>,
                             basic_value_type_t<Vec1>>
Vec1& operator+=(Vec1& v1, const Vec2& v2) {
    auto it = v1.begin(), end = v1.end();
    auto it1 = v2.begin();
    while (it < end) {
        *it += *it1;
        ++it;   ++it1;
    }
    return v1;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires HaveGeneralType<Vec1,Vec2>
auto operator+(const Vec1& v1, const Vec2& v2)
        -> general_vector_type_t<Vec1,Vec2>
{    
    general_vector_type_t<Vec1,Vec2> res(v1.size());
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begin();
    auto it2 = v2.begin();
    while (it < end) {
        *it = *it1 + *it2;
        ++it;   ++it1;  ++it2;
    }
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires GeneralTypeIsFirstOrBoth<Vec1,Vec2>
auto operator+(Vec1&& v1, const Vec2& v2)
        -> std::remove_const_t<std::remove_reference_t<Vec1>>
requires (!std::is_lvalue_reference_v<decltype(v1)>)
{
    std::remove_const_t<std::remove_reference_t<Vec1>> res = std::move(v1);
    res += v2;
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires GeneralTypeIsSecondOrBoth<Vec1,Vec2>
auto operator+(Vec1&& v1, Vec2&& v2)
        -> std::remove_const_t<std::remove_reference_t<Vec2>>
requires (!std::is_lvalue_reference_v<decltype(v2)>)
&& std::is_lvalue_reference_v<decltype(v1)>
{
    return std::move(v2) + v1;
}
//===================
template <VectorLike Vec, NotVectorLike T>
requires std::convertible_to<T,
                             basic_value_type_t<Vec> >
Vec& operator+=(Vec& v, const T& value) {
    auto it = v.begin(), end = v.end();
    while (it < end) {
        *it += value;
        ++it;
    }
    return v;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<Vec,T>
auto operator+(const Vec& vec, const T& value)
        -> vector_t<general_type_t<Vec,T>,vector_dim_v<Vec>>
{
    vector_t<general_type_t<Vec,T>,vector_dim_v<Vec>> res(vec.size());
    auto it  = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it < end) {
        *it = *it1 + value;
        ++it;   ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires GeneralTypeIsFirstOrBoth<Vec,T>
auto operator+(Vec&& v, const T& value)
        -> std::remove_const_t<std::remove_reference_t<Vec>>
requires (!std::is_lvalue_reference_v<decltype(v)>)
{
    std::remove_const_t<std::remove_reference_t<Vec>> res = std::move(v);
    res += value;
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec>
auto operator+(const T& value, Vec&& v)
        -> vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>>
{
    return std::forward<Vec>(v) + value;
}


template <std::input_iterator It1, std::input_iterator It2, VectorLikeIterator It3>
requires ArithmeticIterators<It1,It2,It3>
void sum(It1 first, It2 second, It3 res_first, It3 res_last) {
    while (res_first < res_last) {
        *res_first = *first + *second;
        ++res_first;
        ++first;
        ++second;
    }
}

/*Sum of first and second operands, results is written in range [res_fist,res_last) */
template <std::input_iterator It1, std::input_iterator It2, VectorLikeIteratorOn It3>
requires ArithmeticIterators<It1,It2,It3>
void sum(It1 first, It2 second, It3 res_first, It3 res_last) {
    while (res_first < res_last) {
        sum(first->begin(), second->begin(), res_first->begin(), res_first->end());
        ++res_first;
        ++first;
        ++second;
    }
}

/*==================================
// operator-, -=
==================================*/


template <VectorLike Vec1, VectorsSameDim<Vec1> Vec2>
requires std::convertible_to<basic_value_type_t<Vec2>,
                             basic_value_type_t<Vec1>>
Vec1& operator-=(Vec1& v1, const Vec2& v2) {
    auto it = v1.begin(), end = v1.end();
    auto it1 = v2.begin();
    while (it < end) {
        *it -= *it1;
        ++it;   ++it1;
    }
    return v1;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires HaveGeneralType<Vec1,Vec2>
auto operator-(const Vec1& v1, const Vec2& v2)
        -> general_vector_type_t<Vec1,Vec2>
{
    general_vector_type_t<Vec1,Vec2> res(v1.size());
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begin();
    auto it2 = v2.begin();
    while (it < end) {
        *it = *it1 - *it2;
        ++it;   ++it1;  ++it2;
    }
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires GeneralTypeIsFirstOrBoth<Vec1,Vec2>
auto operator-(Vec1&& v1, const Vec2& v2)
        -> std::remove_const_t<std::remove_reference_t<Vec1>>
requires (!std::is_lvalue_reference_v<decltype(v1)>)
{
    std::remove_const_t<std::remove_reference_t<Vec1>> res = std::move(v1);
    res -= v2;
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires GeneralTypeIsSecondOrBoth<Vec1,Vec2>
auto operator-(Vec1&& v1, Vec2&& v2)
        -> std::remove_const_t<std::remove_reference_t<Vec2>>
requires (!std::is_lvalue_reference_v<decltype(v2)>)
&& std::is_lvalue_reference_v<decltype(v1)>
{
    std::remove_const_t<std::remove_reference_t<Vec2>> res = std::move(v2);
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begin();
    while (it < end) {
        *it = *it1 - *it;
        ++it;   ++it1;
    }
    return res;
}

//===================

template <VectorLike Vec, NotVectorLike T>
requires std::convertible_to<T,
                             basic_value_type_t<Vec> >
Vec& operator-=(Vec& v, const T& value) {
    auto it = v.begin(), end = v.end();
    while (it < end) {
        *it -= value;
        ++it;
    }
    return v;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<Vec,T>
auto operator-(const Vec& vec, const T& value)
        -> vector_t<general_type_t<Vec,T>,vector_dim_v<Vec>>
{
    vector_t<general_type_t<Vec,T>,vector_dim_v<Vec>> res(vec.size());
    auto it  = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it < end) {
        *it = *it1 - value;
        ++it;   ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<Vec,T>
auto operator-(const T& value, const Vec& vec)
        -> vector_t<general_type_t<Vec,T>,vector_dim_v<Vec>>
{
    vector_t<general_type_t<Vec,T>,vector_dim_v<Vec>> res(vec.size());
    auto it  = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it < end) {
        *it = value - *it1;
        ++it;   ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires GeneralTypeIsFirstOrBoth<Vec,T>
auto operator-(Vec&& v, const T& value)
        -> std::remove_const_t<std::remove_reference_t<Vec>>
requires (!std::is_lvalue_reference_v<decltype(v)>)
{
    std::remove_const_t<std::remove_reference_t<Vec>> res = std::move(v);
    res -= value;
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires GeneralTypeIsSecondOrBoth<T,Vec>
auto operator-(const T& value, Vec&& v)
        -> vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>>
requires (!std::is_lvalue_reference_v<decltype(v)>)
{
    std::remove_const_t<std::remove_reference_t<Vec>> res = std::move(v);
    auto it = res.begin(), end = res.end();
    while (it < end) {
        *it = value - *it;
        ++it;
    }
    return res;
}


template <std::input_iterator It1, std::input_iterator It2, VectorLikeIterator It3>
requires ArithmeticIterators<It1,It2,It3>
void sub(It1 first, It2 second, It3 res_first, It3 res_last) {
    while (res_first < res_last) {
        *res_first = *first - *second;
        ++res_first;
        ++first;
        ++second;
    }
}

/*Sum of first and second operands, results is written in range [res_fist,res_last) */
template <std::input_iterator It1, std::input_iterator It2, VectorLikeIteratorOn It3>
requires ArithmeticIterators<It1,It2,It3>
void sub(It1 first, It2 second, It3 res_first, It3 res_last) {
    while (res_first < res_last) {
        sub(first->begin(), second->begin(), res_first->begin(), res_first->end());
        ++res_first;
        ++first;
        ++second;
    }
}

/*==================================
// operator*, *=
==================================*/

template <VectorLike Vec1, VectorsSameDim<Vec1> Vec2>
requires std::convertible_to<basic_value_type_t<Vec2>,
                             basic_value_type_t<Vec1>>
Vec1& operator*=(Vec1& v1, const Vec2& v2) {
    auto it = v1.begin(), end = v1.end();
    auto it1 = v2.begin();
    while (it < end) {
        *it *= *it1;
        ++it;   ++it1;
    }
    return v1;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires HaveGeneralType<Vec1,Vec2>
auto operator*(const Vec1& v1, const Vec2& v2)
        -> general_vector_type_t<Vec1,Vec2>
{    
    general_vector_type_t<Vec1,Vec2> res(v1.size());
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begin();
    auto it2 = v2.begin();
    while (it < end) {
        *it = *it1 * *it2;
        ++it;   ++it1;  ++it2;
    }
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires GeneralTypeIsFirstOrBoth<Vec1,Vec2>
auto operator*(Vec1&& v1, const Vec2& v2)
        -> std::remove_const_t<std::remove_reference_t<Vec1>>
requires (!std::is_lvalue_reference_v<decltype(v1)>)
{
    std::remove_const_t<std::remove_reference_t<Vec1>> res = std::move(v1);
    res *= v2;
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires GeneralTypeIsSecondOrBoth<Vec1,Vec2>
auto operator*(Vec1&& v1, Vec2&& v2)
        -> std::remove_const_t<std::remove_reference_t<Vec2>>
requires (!std::is_lvalue_reference_v<decltype(v2)>)
&& std::is_lvalue_reference_v<decltype(v1)>
{
    return std::move(v2) * v1;
}

//===================

template <VectorLike Vec, NotVectorLike T>
requires std::convertible_to<T,
                             basic_value_type_t<Vec> >
Vec& operator*=(Vec& v, const T& value) {
    auto it = v.begin(), end = v.end();
    while (it < end) {
        *it *= value;
        ++it;
    }
    return v;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<Vec,T>
auto operator*(const Vec& vec, const T& value)
        -> vector_t<general_type_t<Vec,T>,vector_dim_v<Vec>>
{
    vector_t<general_type_t<Vec,T>,vector_dim_v<Vec>> res(vec.size());
    auto it  = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it < end) {
        *it = *it1 * value;
        ++it;   ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires GeneralTypeIsFirstOrBoth<Vec,T>
auto operator*(Vec&& v, const T& value)
        -> std::remove_const_t<std::remove_reference_t<Vec>>
requires (!std::is_lvalue_reference_v<decltype(v)>)
{
    std::remove_const_t<std::remove_reference_t<Vec>> res = std::move(v);
    res *= value;
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec>
auto operator*(const T& value, Vec&& v)
        -> vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>>
{
    return std::forward<Vec>(v) * value;
}


template <std::input_iterator It1, std::input_iterator It2, VectorLikeIterator It3>
requires ArithmeticIterators<It1,It2,It3>
void mul(It1 first, It2 second, It3 res_first, It3 res_last) {
    while (res_first < res_last) {
        *res_first = (*first) * (*second);
        ++res_first;
        ++first;
        ++second;
    }
}

/*Sum of first and second operands, results is written in range [res_fist,res_last) */
template <std::input_iterator It1, std::input_iterator It2, VectorLikeIteratorOn It3>
requires ArithmeticIterators<It1,It2,It3>
void mul(It1 first, It2 second, It3 res_first, It3 res_last) {
    while (res_first < res_last) {
        mul(first->begin(), second->begin(), res_first->begin(), res_first->end());
        ++res_first;
        ++first;
        ++second;
    }
}

/*==================================
// operator/, /=
==================================*/


template <VectorLike Vec1, VectorsSameDim<Vec1> Vec2>
requires std::convertible_to<basic_value_type_t<Vec2>,
                             basic_value_type_t<Vec1>>
Vec1& operator/=(Vec1& v1, const Vec2& v2) {
    auto it = v1.begin(), end = v1.end();
    auto it1 = v2.begin();
    while (it < end) {
        *it /= *it1;
        ++it;   ++it1;
    }
    return v1;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires HaveGeneralType<Vec1,Vec2>
auto operator/(const Vec1& v1, const Vec2& v2)
        -> general_vector_type_t<Vec1,Vec2>
{
    general_vector_type_t<Vec1,Vec2> res(v1.size());
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begin();
    auto it2 = v2.begin();
    while (it < end) {
        *it = *it1 / *it2;
        ++it;   ++it1;  ++it2;
    }
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires GeneralTypeIsFirstOrBoth<Vec1,Vec2>
auto operator/(Vec1&& v1, const Vec2& v2)
        -> std::remove_const_t<std::remove_reference_t<Vec1>>
requires (!std::is_lvalue_reference_v<decltype(v1)>)
{
    std::remove_const_t<std::remove_reference_t<Vec1>> res = std::move(v1);
    res /= v2;
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires GeneralTypeIsSecondOrBoth<Vec1,Vec2>
auto operator/(Vec1&& v1, Vec2&& v2)
        -> std::remove_const_t<std::remove_reference_t<Vec2>>
requires (!std::is_lvalue_reference_v<decltype(v2)>)
&& std::is_lvalue_reference_v<decltype(v1)>
{
    std::remove_const_t<std::remove_reference_t<Vec2>> res = std::move(v2);
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begin();
    while (it < end) {
        *it = *it1 / *it;
        ++it;   ++it1;
    }
    return res;
}

//===================

template <VectorLike Vec, NotVectorLike T>
requires std::convertible_to<T,
                             basic_value_type_t<Vec> >
Vec& operator/=(Vec& v, const T& value) {
    auto it = v.begin(), end = v.end();
    while (it < end) {
        *it /= value;
        ++it;
    }
    return v;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<Vec,T>
auto operator/(const Vec& vec, const T& value)
        -> vector_t<general_type_t<Vec,T>,vector_dim_v<Vec>>
{
    vector_t<general_type_t<Vec,T>,vector_dim_v<Vec>> res(vec.size());
    auto it  = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it < end) {
        *it = *it1 / value;
        ++it;   ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<Vec,T>
auto operator/(const T& value, const Vec& vec)
        -> vector_t<general_type_t<Vec,T>,vector_dim_v<Vec>>
{
    vector_t<general_type_t<Vec,T>,vector_dim_v<Vec>> res(vec.size());
    auto it  = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it < end) {
        *it = value / *it1;
        ++it;   ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires GeneralTypeIsFirstOrBoth<Vec,T>
auto operator/(Vec&& v, const T& value)
        -> std::remove_const_t<std::remove_reference_t<Vec>>
requires (!std::is_lvalue_reference_v<decltype(v)>)
{
    std::remove_const_t<std::remove_reference_t<Vec>> res = std::move(v);
    res /= value;
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires GeneralTypeIsSecondOrBoth<T,Vec>
auto operator/(const T& value, Vec&& v)
        -> vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>>
requires (!std::is_lvalue_reference_v<decltype(v)>)
{
    std::remove_const_t<std::remove_reference_t<Vec>> res = std::move(v);
    auto it = res.begin(), end = res.end();
    while (it < end) {
        *it = value / *it;
        ++it;
    }
    return res;
}


template <std::input_iterator It1, std::input_iterator It2, VectorLikeIterator It3>
requires ArithmeticIterators<It1,It2,It3>
void div(It1 first, It2 second, It3 res_first, It3 res_last) {
    while (res_first < res_last) {
        *res_first = *first / *second;
        ++res_first;
        ++first;
        ++second;
    }
}

/*Sum of first and second operands, results is written in range [res_fist,res_last) */
template <std::input_iterator It1, std::input_iterator It2, VectorLikeIteratorOn It3>
requires ArithmeticIterators<It1,It2,It3>
void div(It1 first, It2 second, It3 res_first, It3 res_last) {
    while (res_first < res_last) {
        div(first->begin(), second->begin(), res_first->begin(), res_first->end());
        ++res_first;
        ++first;
        ++second;
    }
}

/*==================================
// Other methods/operators overloadings
==================================*/

// operator<<
template<VectorLike Vec>
requires std::same_as<typename Vec::value_type
                    , typename Vec::basic_value_type>
std::ostream& operator<<(std::ostream& stream, const Vec& v) {
    auto it = v.begin(), end = v.end();
    while (it < end) {
        if (*it >= 0) stream << ' ';
        stream << *it << ' ';
        ++it;
    }
    return stream;
}
// overloading to print vector Dim > 1 by rows
template<VectorLike Vec>
std::ostream& operator<<(std::ostream& stream, const Vec& v) {
    auto it = v.begin(), end = v.end();
    while (it < end) {
        stream << *it;
        if ((*it).size())
            stream << "\n";
        ++it;
    }
    return stream;
}

namespace detail {

// get vector_t<T,Dim> size as std::array<size_t,Dim>
template <VectorLike V>
void size_helper(const V& vec, size_t* sizes) {
    *sizes = vec.size();
    if constexpr (is_vector_like_v<typename V::value_type>) {
        if (vec.size() > 0) {
            size_helper<typename V::value_type>(vec[0], ++sizes);
        }
    }
}

} // namespace detail

template <VectorLike V>
auto size(const V& vec) 
        -> std::array<size_t,vector_dim_v<V> > {
    const size_t dim = vector_dim_v<V>;
    std::array<size_t,dim> sizes{};
    detail::size_helper<V>(vec,sizes.begin());
    return sizes;
}


template <typename T>
void swap(vector<T>& lhs, vector<T>& rhs) {
    std::swap(lhs.v, rhs.v);
}

template <VectorLike V1, VectorLike V2>
requires std::same_as<typename V1::value_type, typename V2::value_type>
&& VectorLike<typename V1::value_type>
&& (!std::same_as<V1,V2>)
void swap(V1& lhs, V2& rhs) {
    auto itl = lhs.begin(), itl_end = lhs.end();
    auto itr = rhs.begin();
    while (itl < itl_end) {
        swap(*itl,*itr);
        ++itl;
        ++itr;
    }
}

/*For vectors*/
template <VectorLikeIterator It, NotVectorLike T>
void fill(It first, It last, const T& value) {
    std::fill(first,last,value);
}

/*For matrices (2d vectors)*/
template <VectorLikeIteratorOn It, NotVectorLike T>
void fill(It first, It last, const T& value) {
    while (first < last) {
        fill(first->begin(), first->end(), value);
        ++first;
    }
}

/* Fill first n1,n2,... elements of VectorLike */
namespace detail {

template <VectorLike V, typename T>
requires std::is_convertible_v<T,typename V::value_type>
void fill_helper(V& v, const T& value, size_t sz) {
    fill(v.begin(),v.begin()+sz,value);
}

template <VectorLike V, typename T, typename... Sizes>
void fill_helper(V& v, const T& value, size_t sz, Sizes... sizes) {
    auto it = v.begin();
    auto end = it + sz;
    while (it < end) {
        fill_helper(*it,value,sizes...);
        ++it;
    }
}


} // namespace detail


template <VectorLike V, typename T, typename... Sizes>
requires (std::is_convertible_v<Sizes,size_t> && ...)
void fill(V& v, const T& value, Sizes... sizes) {
    detail::fill_helper(v,value,sizes...);
}


template <VectorLike V, typename T>
requires std::is_convertible_v<T,typename V::value_type>
void fill(V& v, const T& value) {
    fill(v.begin(),v.end(),value);
}

template <VectorLike V, typename T>
void fill(V& v, const T& value) {
    auto it = v.begin();
    auto end = v.end();
    while (it < end) {
        fill(*it,value);
        ++it;
    }
}


} // namespace math 