#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <concepts>
#include <exception>
#include <type_traits>
#include <concepts>


namespace math {

/*==================================
// class vector declaration
==================================*/
template <typename T>
class vector;

template <typename It1, typename It2>
requires std::same_as<typename It1::value_type
                    , typename It2::value_type>
class Slice;

/*==================================
// vector meta-functions and concepts
==================================*/

template <typename T, size_t Dim = 1>
struct multidim_vector {
    using type = vector< typename multidim_vector<T, Dim-1>::type >;
};
template<typename T>
struct multidim_vector<T,0> {
    using type = T;
};
template <typename T, size_t Dim = 1>
using vector_t = typename multidim_vector<T,Dim>::type;


template <typename T>
concept VectorLike = requires() {
    std::remove_reference_t<T>::is_vector_like;
};

template <typename T>
struct is_vector_like {
    static const bool value = false;
};
template <VectorLike T>
struct is_vector_like<T> {
    static const bool value = true;
};
template <typename T>
static const bool is_vector_like_v = is_vector_like<T>::value;

template <typename T>
concept NotVectorLike = (!VectorLike<T>);

template <typename vec>
struct vector_dim {
    static const size_t value = (vector_dim<typename std::remove_reference_t<vec>::value_type>::value + 1);
};
template <NotVectorLike T>
struct vector_dim<T> {
    static const size_t value = 0;
};
template <VectorLike T>
static const size_t vector_dim_v = vector_dim<T>::value;

template <typename T>
struct basic_value_type_helper {
    using type = T;
};

template <VectorLike vec>
struct basic_value_type_helper<vec> {
    using type = typename basic_value_type_helper<typename vec::value_type>::type;
};



template <typename T>
using basic_value_type_t = typename basic_value_type_helper<
                                std::remove_const_t<std::remove_reference_t<T>>>::type;


template <typename Vec1, typename Vec2>
concept VectorsSameDim = VectorLike<Vec1>
                    &&   VectorLike<Vec2>
                    &&   vector_dim_v<Vec1> == vector_dim_v<Vec2>;

template <typename T>
concept Vector = VectorLike<T>
              && NotVectorLike<typename T::value_type>;

template <typename T>
concept Matrix = VectorLike<T>
              && Vector<typename T::value_type>;


template <typename T, typename U>
struct general_type {
    using type = std::remove_reference_t<decltype(true
        ? std::declval<basic_value_type_t<T>>()
        : std::declval<basic_value_type_t<U>>())>;
};
template <typename T, typename U>
using general_type_t = typename general_type<T,U>::type;


template <typename T, typename U>
concept HaveGeneralType = requires() {
    typename general_type_t<T,U>;    
};

template <typename First, typename Second>
concept GeneralTypeIsOther =
    HaveGeneralType<First,Second>
    && (!std::same_as<general_type_t<First,Second>,basic_value_type_t<First>>)
    && (!std::same_as<general_type_t<First,Second>,basic_value_type_t<Second>>);

template <typename First, typename Second>
concept GeneralTypeIsFirst =
    HaveGeneralType<First,Second>
    &&   std::same_as<general_type_t<First,Second>,basic_value_type_t<First>>
    && (!std::same_as<general_type_t<First,Second>,basic_value_type_t<Second>>);

template <typename First, typename Second>
concept GeneralTypeIsSecond =
    HaveGeneralType<First,Second>
    && (!std::same_as<general_type_t<First,Second>,basic_value_type_t<First>>)
    &&   std::same_as<general_type_t<First,Second>,basic_value_type_t<Second>>;

template <typename First, typename Second>
concept GeneralTypeIsFirstOrBoth =
    GeneralTypeIsFirst<First,Second>
    || std::same_as<basic_value_type_t<First>,
                    basic_value_type_t<Second>>;

template <typename First, typename Second>
concept GeneralTypeIsSecondOrBoth =
    GeneralTypeIsSecond<First,Second>
    || std::same_as<basic_value_type_t<First>,
                    basic_value_type_t<Second>>;

template <typename First, typename Second>
concept GeneralTypeIsFirstOrOther =
    GeneralTypeIsFirst<First,Second>
    || GeneralTypeIsOther<First,Second>;

template <typename First, typename Second>
concept GeneralTypeIsSecondOrOther =
    GeneralTypeIsSecond<First,Second>
    || GeneralTypeIsOther<First,Second>;


template <typename Vec1, typename Vec2>
concept ArithmeticVectorsLike = HaveGeneralType<Vec1,Vec2>
&& VectorsSameDim<Vec1,Vec2>;


template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
struct general_vector_type {
    using type = vector_t<general_type_t<Vec1,Vec2>
                        , vector_dim_v<Vec1>>;
};
template <typename Vec1, typename Vec2>
using general_vector_type_t = typename general_vector_type<Vec1,Vec2>::type;


/*==================================
// class vector definition
==================================*/

template <typename T>
class vector {
    std::vector<T> v;
public:
    using value_type = T;
    using reference = T&;
    static const bool is_vector_like = true;
    using basic_value_type = basic_value_type_t<T>;
    
    
    using iterator       = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    
    constexpr iterator       begin();
    constexpr const_iterator begin() const;
    constexpr const_iterator cbegin() const;
    constexpr iterator       end();
    constexpr const_iterator end()   const;
    constexpr const_iterator cend()   const;
    
    
    vector()                             : v()        {}//{ std::cout << "default\n";}
    vector(std::initializer_list<T> init): v(init)    {}//{ std::cout << "init\n"; }
    vector(size_t size)                  : v(size)    {}//{ std::cout << "size_t\n";;}
    vector(const vector<T>& other)       : v(other.v) {}//{ std::cout << "copy\n";; }

    vector(vector&& other) noexcept      : v(std::move(other.v)) {}//{ std::cout << "move\n";; }
    vector(const std::vector<T> other)   : v(other) {}//{std::cout << "std::vector\n";; }
    
    vector(size_t size, const T& value)   : v(size,value) {}//{ std::cout << "size_t, value\n";}

    template <std::input_iterator InputIt>
    vector(InputIt first, InputIt last) : v(first,last) {}

    vector(const Slice<iterator,iterator>& sl);
    vector(const Slice<const_iterator,const_iterator>& sl);

    ~vector() {}//{ std::cout << "~vector [" << *this << "]"<< std::endl; };

    vector<T>& operator=(const vector<T>&) &;
    vector<T>& operator=(vector<T>&&) & ;
    vector<T>& operator=(const Slice<iterator,iterator>&) &;
    vector<T>& operator=(const Slice<const_iterator,const_iterator>&) &;

    template <ArithmeticVectorsLike<vector<T>> V2>
    vector<T>& operator=(const V2&) &;
    
    T& operator[](size_t) &;
    const T& operator[](size_t) const &;
    T&& operator[](size_t) &&;
    const T&& operator[](size_t) const &&;
    T& last() &;
    const T& last() const &;
    T&& last() &&;
    const T&& last() const &&;

    // auto&& operator[](this auto&& self, size_t pos) {
    //     return std::forward_like<decltype(self)>(self.v[pos])
    // }
    T& at(size_t);
    const T& at(size_t) const;


    
    size_t size() const;
    
    void push_back(const T&);
    void push_back(T&&);

    template <typename... Args>
    void emplace_back(Args&&...); // universal reference

    


};

template <typename It1, typename It2>
vector(const Slice<It1,It2>&) -> vector<typename It1::value_type>;

template <typename T>
vector<T>::vector(const Slice<iterator, iterator> &sl) : v(sl.size()) {
    auto it_ = begin();
    auto it1 = sl.begin(), end1 = sl.end();
    while (it1 < end1) {
        *it_ = *it1;
        ++it_;   ++it1;
    }
}
template <typename T>
vector<T>::vector(const Slice<const_iterator, const_iterator> &sl) : v(sl.size()) {
    auto it_ = begin();
    auto it1 = sl.begin(), end1 = sl.end();
    while (it1 < end1) {
        *it_ = *it1;
        ++it_;   ++it1;
    }
}


/*==================================
// class vector methods
==================================*/

/*==================================
// class std::vector methods
==================================*/

template <typename T>
size_t vector<T>::size() const {
    return v.size();
}

template <typename T>
void vector<T>::push_back(const T& value) {
    v.push_back(value);
}
template <typename T>
void vector<T>::push_back(T&& value) {
    v.push_back(std::move(value));
}

template <typename T>
template <typename... Args>
void vector<T>::emplace_back(Args&&... args) {
    v.emplace_back(std::forward<Args>(args)...);
}



// operator=

template <typename T>
vector<T>& vector<T>::operator=(const vector<T>& other) & {
    v = other.v;
    return *this;
}

template <typename T>
vector<T>& vector<T>::operator=(vector<T>&& other) & {
    v = std::move(other.v);
    return *this;
}

template <typename T>
template <ArithmeticVectorsLike<vector<T>> V2>
vector<T>& vector<T>::operator=(const V2& other) & {
    auto it = begin(), end_ = end();
    auto it1 = other.begin();
    while (it < end_) {
        *it = *it1;
        ++it;   ++it1;
    }
    return *this;
};

template <typename T>
vector<T>& vector<T>::operator=(const Slice<iterator, iterator>& sl) & {
    v = std::vector<T>(sl.size());
    auto it = begin();
    auto it1 = sl.begin(), end1 = sl.end();
    while (it1 < end1) {
        *it = *it1;
        ++it; ++it1;
    }
    return *this;
}

template <typename T>
vector<T>& vector<T>::operator=(const Slice<const_iterator, const_iterator>& sl) & {
    v = std::vector<T>(sl.size());
    auto it = begin();
    auto it1 = sl.begin(), end1 = sl.end();
    while (it1 < end1) {
        *it = *it1;
        ++it; ++it1;
    }
    return *this;
}

//operator[]

template <typename T>
T& vector<T>::operator[](size_t pos) & {
    return v[pos];
}
template <typename T>
const T& vector<T>::operator[](size_t pos) const & {
    return v[pos];
}
template <typename T>
T&& vector<T>::operator[](size_t pos) && {
    return std::move(v[pos]);
}
template <typename T>
const T&& vector<T>::operator[](size_t pos) const && {
    return std::move(v[pos]);
}

// template <typename T>
// auto&& vector<T>::operator[](this auto&& self, size_t pos) {
//     return std::forward_like<decltype(self)>(v[pos]);
// }

template <typename T>
T& vector<T>::last() & {
    return v[size()-1];
}
template <typename T>
const T& vector<T>::last() const & {
    return v[size()-1];
}
template <typename T>
T&& vector<T>::last() && {
    return std::move(v[size()-1]);
}
template <typename T>
const T&& vector<T>::last() const && {
    return std::move(v[size()-1]);
}

template <typename T>
T& vector<T>::at(size_t pos) {
    return v.at(pos);
}

template <typename T>
const T& vector<T>::at(size_t pos) const {
    return v.at(pos);
}

// get iterator methods

template<typename T>
constexpr typename vector<T>::iterator vector<T>::begin() {
    return v.begin();
} 

template<typename T>
constexpr typename vector<T>::const_iterator vector<T>::begin() const {
    return v.begin();
}

template<typename T>
constexpr typename vector<T>::const_iterator vector<T>::cbegin() const {
    return v.cbegin();
} 

template<typename T>
constexpr typename vector<T>::iterator vector<T>::end() {
    return v.end();
} 

template<typename T>
constexpr typename vector<T>::const_iterator vector<T>::end() const {
    return v.end();
}

template<typename T>
constexpr typename vector<T>::const_iterator vector<T>::cend() const {
    return v.cend();
}



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

// get vector_t<T,Dim> size as std::array<size_t,Dim>
template <typename T>
void size_helper(const vector<T>& vec, size_t* sizes) {
    *sizes = vec.size();
    if constexpr (is_vector_like_v<T>) {
        if (vec.size() > 0) {
            size_helper<typename T::value_type>(vec[0], ++sizes);
        }
    }
}

template <typename T>
auto size(const vector<T>& vec) 
        -> std::array<size_t,vector_dim_v<vector<T>> > {
    const size_t dim = vector_dim_v<vector<T>>;
    std::array<size_t,dim> sizes{};
    size_helper<T>(vec,sizes.begin());
    return sizes;
}



} // namespace math