#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <concepts>
#include <exception>
#include <type_traits>

namespace math {


/*==================================
// class vector declaration
==================================*/
template <typename T>
class vector;

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
    T::is_vector_like;
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
concept NotVectorLike = (!is_vector_like_v<T>);

template <typename vec>
struct vector_dim {
    static const size_t value = (vector_dim<typename vec::value_type>::value + 1);
};
template <typename T>
requires (!is_vector_like_v<T>)
struct vector_dim<T> {
    static const size_t value = 0;
};
template <VectorLike T>
static const size_t vector_dim_v = vector_dim<T>::value;

template <typename vec>
struct vector_basic_value_type_helper {
    using type = typename vector_basic_value_type_helper<typename vec::value_type>::type;
};

template <typename T>
requires (!is_vector_like_v<T>)
struct vector_basic_value_type_helper<T> {
    using type = T;
};

template <typename vec>
requires (is_vector_like_v<vec>)
struct vector_basic_value_type {
    using type = typename vector_basic_value_type_helper<vec>::type;
};

template <typename vec>
using vector_basic_value_type_t = typename vector_basic_value_type<vec>::type;



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
    using type = decltype(true
                    ? std::declval<T>()
                    : std::declval<U>()); 
};

template <VectorLike T, NotVectorLike U>
struct general_type<T,U> {
    using type = std::decay_t<decltype(true
        ? std::declval<typename T::basic_value_type>()
        : std::declval<U>())>;
};
template <NotVectorLike T, VectorLike U>
struct general_type<T,U> {
    using type = std::decay_t<decltype(true
        ? std::declval<T>()
        : std::declval<typename U::basic_value_type>())>;
};
template <VectorLike T, VectorLike U>
struct general_type<T,U> {
    using type = std::decay_t<decltype(true
        ? std::declval<typename T::basic_value_type>()
        : std::declval<typename U::basic_value_type>())>;
};
template <typename T, typename U>
using general_type_t = typename general_type<T,U>::type;

template <typename T, typename U>
concept HaveGeneralType = requires() {
    typename general_type_t<T,U>;    
};

template <typename Vec1, typename Vec2>
concept ArithmeticVectorsLike = HaveGeneralType<Vec1,Vec2>
&& VectorsSameDim<Vec1,Vec2>;


template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
struct general_arithmetic_vector_like_type {
    using type = vector_t<general_type_t<Vec1,Vec2>
                          , vector_dim_v<Vec1>>;
};
template <typename Vec1, typename Vec2>
using general_arithmetic_vector_like_type_t = typename general_arithmetic_vector_like_type<Vec1,Vec2>::type;

/*==================================
// class vector definition
==================================*/

template <typename T>
class vector {
    std::vector<T> v;
public:
    using value_type = T;
    using reference = value_type&;
    static const bool is_vector_like = true;
    using basic_value_type = typename vector_basic_value_type_helper<T>::type;
    
    
    using iterator       = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    
    iterator       begin();
    const_iterator begin() const;
    iterator       end();
    const_iterator end()   const;
    
    
    vector()                             : v()        {}//{ std::cout << "default\n";}
    vector(std::initializer_list<T> init): v(init)    {}//{ std::cout << "init\n"; }
    vector(size_t size)                  : v(size)    {}//{ std::cout << "size_t\n";;}
    vector(const vector<T>& other)       : v(other.v) {}//{ std::cout << "copy\n";; }

    vector(vector&& other) noexcept      : v(std::move(other.v)) {}//{ std::cout << "move\n";; }
    vector(const std::vector<T> other)   : v(other) {}//{std::cout << "std::vector\n";; }
    
    vector(size_t size, const T& value)   : v(size,value) {}//{ std::cout << "size_t, value\n";}

    ~vector() {}//{ std::cout << "~vector [" << *this << "]"<< std::endl; };

    vector<T>& operator=(const vector<T>&) &;
    vector<T>& operator=(vector<T>&&) &;
    
    T& operator[](size_t) &;
    const T& operator[](size_t) const &;
    T&& operator[](size_t) &&;
    const T&& operator[](size_t) const &&;

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


/*==================================
// class vector methods
==================================*/

/*==================================
// class std::vector methods
==================================*/

template <typename T>
size_t vector<T>::size() const {
    return this->v.size();
}

template <typename T>
void vector<T>::push_back(const T& value) {
    this->v.push_back(value);
}
template <typename T>
void vector<T>::push_back(T&& value) {
    this->v.push_back(std::move(value));
}

template <typename T>
template <typename... Args>
void vector<T>::emplace_back(Args&&... args) {
    this->v.emplace_back(std::forward<Args>(args)...);
}

// operator=

template <typename T>
vector<T>& vector<T>::operator=(const vector<T>& other) & {
    this->v = other.v;
    return *this;
}

template <typename T>
vector<T>& vector<T>::operator=(vector<T>&& other) & {
    this->v = std::move(other.v);
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
T& vector<T>::at(size_t pos) {
    return this->v.at(pos);
}

template <typename T>
const T& vector<T>::at(size_t pos) const {
    return this->v.at(pos);
}

// get iterator methods

template<typename T>
typename vector<T>::iterator vector<T>::begin() {
    return this->v.begin();
} 

template<typename T>
typename vector<T>::const_iterator vector<T>::begin() const {
    return this->v.begin();
} 

template<typename T>
typename vector<T>::iterator vector<T>::end() {
    return this->v.end();
} 

template<typename T>
typename vector<T>::const_iterator vector<T>::end() const {
    return this->v.end();
} 



/*==================================
// vector mathematical operators overloadings
==================================*/

/*==================================
// operator+, +=
==================================*/

template <VectorLike Vec1, VectorsSameDim<Vec1> Vec2>
requires std::is_convertible_v<typename Vec2::basic_value_type,
                               typename Vec1::basic_value_type>
Vec1& operator+=(Vec1& v1, const Vec2& v2) {
    auto it = v1.begin(), end = v1.end();
    auto it1 = v2.begin();
    while (it != end) {
        *it += *it1;
        ++it;   ++it1;
    }
    return v1;
}

template <VectorLike Vec, NotVectorLike T>
requires std::is_convertible_v<T,
                               typename Vec::basic_value_type>
Vec& operator+=(Vec& v, const T& value) {
    auto it = v.begin(), end = v.end();
    while (it != end) {
        *it += value;
        ++it;
    }
    return v;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
auto operator+(const Vec1& v1, const Vec2& v2)
        -> general_arithmetic_vector_like_type_t<Vec1,Vec2> {
    
    general_arithmetic_vector_like_type_t<Vec1,Vec2> res = zeros<general_type_t<Vec1,Vec2>>(v1);
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begin();
    auto it2 = v2.begin();
    while (it != end) {
        *it = *it1 + *it2;
        ++it;   ++it1;  ++it2;
    }
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec1>
Vec1 operator+(const Vec1& v1, const Vec2& v2) {
    Vec1 res = v1;
    res += v2;
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec1>
Vec1 operator+(Vec1&& v1, const Vec2& v2) {
    Vec1 res = std::move(v1);
    res += v2;
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec2>
Vec2 operator+(const Vec1& v1, const Vec2& v2) {
    return v2 + v1;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec2>
Vec2 operator+(const Vec1& v1, Vec2&& v2) {
    return std::move(v2) + v1;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
auto operator+(const Vec& vec, const T& value)
        -> vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> {
    vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> res = zeros<general_type_t<T,Vec>>(vec);
    auto it  = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it != end) {
        *it = (*it1) + value;
        ++it;
        ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
auto operator+(const T& value, const Vec& vec)
        -> vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> {
    return vec + value;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
&& std::is_same_v<general_type_t<Vec,T>,typename Vec::basic_value_type>
Vec operator+(const Vec& v, const T& value) {
    Vec res = v;
    res += value;
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
&& std::is_same_v<general_type_t<Vec,T>,typename Vec::basic_value_type>
Vec operator+(Vec&& v, const T& value) {
    Vec res = std::move(v);
    res += value;
    return res;
}


/*==================================
// operator-, -=
==================================*/


template <VectorLike Vec1, VectorsSameDim<Vec1> Vec2>
requires std::is_convertible_v<typename Vec2::basic_value_type,
                               typename Vec1::basic_value_type>
Vec1& operator-=(Vec1& v1, const Vec2& v2) {
    auto it = v1.begin(), end = v1.end();
    auto it1 = v2.begin();
    while (it != end) {
        *it -= *it1;
        ++it;   ++it1;
    }
    return v1;
}

template <VectorLike Vec, NotVectorLike T>
requires std::is_convertible_v<T,
                               typename Vec::basic_value_type>
Vec& operator-=(Vec& v, const T& value) {
    auto it = v.begin(), end = v.end();
    while (it != end) {
        *it -= value;
        ++it;
    }
    return v;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
auto operator-(const Vec1& v1, const Vec2& v2)
        -> general_arithmetic_vector_like_type_t<Vec1,Vec2> {
    general_arithmetic_vector_like_type_t<Vec1,Vec2> res = zeros<general_type_t<Vec1,Vec2>>(v1);
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begin();
    auto it2 = v2.begin();
    while (it != end) {
        *it = *it1 - *it2;
        ++it;   ++it1;  ++it2;
    }
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec1>
Vec1 operator-(const Vec1& v1, const Vec2& v2) {
    Vec1 res = v1;
    res -= v2;
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec1>
Vec1 operator-(Vec1&& v1, const Vec2& v2) {
    Vec1 res = std::move(v1);
    res -= v2;
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec2>
Vec2 operator-(const Vec1& v1, const Vec2& v2) {
    Vec2 res = v2;
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begin();
    while (it != end) {
        *it = *it1 - *it;
        ++it;   ++it1;
    }
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec2>
Vec2 operator-(const Vec1& v1, Vec2&& v2) {
    Vec2 res = std::move(v2);
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begin();
    while (it != end) {
        *it = *it1 - *it;
        ++it;   ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
auto operator-(const Vec& vec, const T& value)
        -> vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> {
    vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> res = zeros<general_type_t<T,Vec>>(vec);
    auto it  = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it != end) {
        *it = (*it1) - value;
        ++it;
        ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
auto operator-(const T& value, const Vec& vec)
        -> vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> {
    vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> res = zeros<general_type_t<T,Vec>>(vec);
    auto it  = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it != end) {
        *it =  value - (*it1);
        ++it;
        ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
&& std::is_same_v<general_type_t<Vec,T>,typename Vec::basic_value_type>
Vec operator-(const Vec& v, const T& value) {
    Vec res = v;
    res -= value;
    return res;
}

/*==================================
// operator*, *=
==================================*/


template <VectorLike Vec1, VectorsSameDim<Vec1> Vec2>
requires std::is_convertible_v<typename Vec2::basic_value_type,
                               typename Vec1::basic_value_type>
Vec1& operator*=(Vec1& v1, const Vec2& v2) {
    auto it = v1.begin(), end = v1.end();
    auto it1 = v2.begin();
    while (it != end) {
        *it *= *it1;
        ++it;   ++it1;
    }
    return v1;
}

template <VectorLike Vec, NotVectorLike T>
requires std::is_convertible_v<T,
                               typename Vec::basic_value_type>
Vec& operator*=(Vec& v, const T& value) {
    auto it = v.begin(), end = v.end();
    while (it != end) {
        *it *= value;
        ++it;
    }
    return v;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
auto operator*(const Vec1& v1, const Vec2& v2)
        -> general_arithmetic_vector_like_type_t<Vec1,Vec2> {
    general_arithmetic_vector_like_type_t<Vec1,Vec2> res = zeros<general_type_t<Vec1,Vec2>>(v1);
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begin();
    auto it2 = v2.begin();
    while (it != end) {
        *it = *it1 * *it2;
        ++it;   ++it1;  ++it2;
    }
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec1>
Vec1 operator*(const Vec1& v1, const Vec2& v2) {
    Vec1 res = v1;
    res *= v2;
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec1>
Vec1 operator*(Vec1&& v1, const Vec2& v2) {
    Vec1 res = std::move(v1);
    res *= v2;
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec2>
Vec2 operator*(const Vec1& v1, const Vec2& v2) {
    return v2 * v1;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec2>
Vec2 operator*(const Vec1& v1, Vec2&& v2) {
    return std::move(v2) * v1;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
auto operator*(const Vec& vec, const T& value)
        -> vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> {
    vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> res = zeros<general_type_t<T,Vec>>(vec);
    auto it  = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it != end) {
        *it = (*it1) * value;
        ++it;
        ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
auto operator*(const T& value, const Vec& vec)
        -> vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> {
    return vec * value;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
&& std::is_same_v<general_type_t<Vec,T>,typename Vec::basic_value_type>
Vec operator*(const Vec& v, const T& value) {
    Vec res = v;
    res *= value;
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
&& std::is_same_v<general_type_t<Vec,T>,typename Vec::basic_value_type>
Vec operator*(Vec&& v, const T& value) {
    Vec res = std::move(v);
    res *= value;
    return res;
}

/*==================================
// operator/, +/=
==================================*/

template <VectorLike Vec1, VectorsSameDim<Vec1> Vec2>
requires std::is_convertible_v<typename Vec2::basic_value_type,
                               typename Vec1::basic_value_type>
Vec1& operator/=(Vec1& v1, const Vec2& v2) {
    auto it = v1.begin(), end = v1.end();
    auto it1 = v2.begin();
    while (it != end) {
        *it /= *it1;
        ++it;   ++it1;
    }
    return v1;
}

template <VectorLike Vec, NotVectorLike T>
requires std::is_convertible_v<T,
                               typename Vec::basic_value_type>
Vec& operator/=(Vec& v, const T& value) {
    auto it = v.begin(), end = v.end();
    while (it != end) {
        *it /= value;
        ++it;
    }
    return v;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
auto operator/(const Vec1& v1, const Vec2& v2)
        -> general_arithmetic_vector_like_type_t<Vec1,Vec2> {
    general_arithmetic_vector_like_type_t<Vec1,Vec2> res = zeros<general_type_t<Vec1,Vec2>>(v1);
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begin();
    auto it2 = v2.begin();
    while (it != end) {
        *it = *it1 / *it2;
        ++it;   ++it1;  ++it2;
    }
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec1>
Vec1 operator/(const Vec1& v1, const Vec2& v2) {
    Vec1 res = v1;
    res /= v2;
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec1>
Vec1 operator/(Vec1&& v1, const Vec2& v2) {
    Vec1 res = std::move(v1);
    res /= v2;
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec2>
Vec2 operator/(const Vec1& v1, const Vec2& v2) {
    Vec2 res = v2;
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begini();
    while (it != end) {
        *it = *it1 / *it;
        ++it;   ++it1;
    }
    return res;
}

template <VectorLike Vec1, ArithmeticVectorsLike<Vec1> Vec2>
requires std::is_same_v<general_arithmetic_vector_like_type_t<Vec1,Vec2>,Vec2>
Vec2 operator/(const Vec1& v1, Vec2&& v2) {
    Vec2 res = std::move(v2);
    auto it = res.begin(), end = res.end();
    auto it1 = v1.begini();
    while (it != end) {
        *it = *it1 / *it;
        ++it;   ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
auto operator/(const Vec& vec, const T& value)
        -> vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> {
    vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> res = zeros<general_type_t<T,Vec>>(vec);
    auto it  = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it != end) {
        *it = (*it1) / value;
        ++it;
        ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
auto operator/(const T& value, const Vec& vec)
        -> vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> {
    vector_t<general_type_t<T,Vec>,vector_dim_v<Vec>> res = zeros<general_type_t<T,Vec>>(vec);
    auto it  = res.begin(), end = res.end();
    auto it1 = vec.begin();
    while (it != end) {
        *it =  value / (*it1);
        ++it;
        ++it1;
    }
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
&& std::is_same_v<general_type_t<Vec,T>,typename Vec::basic_value_type>
Vec operator/(const Vec& v, const T& value) {
    Vec res = v;
    res /= value;
    return res;
}

template <VectorLike Vec, NotVectorLike T>
requires HaveGeneralType<T,Vec> 
&& std::is_same_v<general_type_t<Vec,T>,typename Vec::basic_value_type>
Vec operator/(Vec&& v, const T& value) {
    Vec res = std::move(v);
    res /= value;
    return res;
}

/*==================================
// Other methods/operators overloadings
==================================*/

// operator<<
template<typename T>
std::ostream& operator<<(std::ostream& stream, const vector<T>& v) {
    for (const T& x: v) {
        stream << x << " ";
    }
    return stream;
}
// overloading to print vector Dim > 1 by rows
template<typename T>
std::ostream& operator<<(std::ostream& stream, const vector<vector<T>>& v) {
    for (const vector<T>& x: v) {
        stream << x << "\n";
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









template <typename T>
struct Debug {
    Debug() = delete;
};

} // namespace math