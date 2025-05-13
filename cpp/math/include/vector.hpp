#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <concepts>
// #include <exception>
#include <type_traits>
#include <concepts>

#include "vectorlike_type_traits.hpp"

namespace math {

/*==================================
// class vector declaration
==================================*/
template <typename T>
class vector;

template <typename It1, typename It2>
requires ::std::same_as<typename It1::value_type
                      , typename It2::value_type>
class Slice;


/*=================================

===================================*/

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

    explicit vector(const Slice<iterator,iterator>& sl);
    explicit vector(const Slice<const_iterator,const_iterator>& sl);

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

    void resize(size_t count);
    void resize(size_t count, const T& value);

    iterator erase(const_iterator pos);
    iterator erase (const_iterator first, const_iterator last);

    template <typename U>
    friend void swap(vector<U>&, vector<U>&);

    T* data();
    const T* data() const;
};

template <typename It1, typename It2>
vector(const Slice<It1,It2>&) -> vector<typename It1::value_type>;

template <typename T>
vector<T>::vector(const Slice<iterator, iterator> &sl)
        : vector<T>(sl.begin(),sl.end()) {
}
template <typename T>
vector<T>::vector(const Slice<const_iterator, const_iterator> &sl)
        : vector<T>(sl.begin(),sl.end()) {
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

template <typename T>
void math::vector<T>::resize(size_t count) {
    v.resize(count);
}

template <typename T>
void math::vector<T>::resize(size_t count, const T& value) {
    v.resize(count,value);
}

template <typename T>
vector<T>::iterator vector<T>::erase(typename vector<T>::const_iterator pos) {
    return v.erase(pos);
}
template <typename T>
vector<T>::iterator vector<T>::erase(typename vector<T>::const_iterator first, typename vector<T>::const_iterator last) {
    return v.erase(first,last);
}


template <typename T>
T* vector<T>::data() {
    return v.data();
}


template <typename T>
const T* vector<T>::data() const {
    return v.data();
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







} // namespace math