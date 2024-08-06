#include <iostream>
#include <vector>
#include <algorithm>
#include <concepts>
#include <variant>

namespace math {

template <typename T>
class vector {
    
public:
    using value_type = T;
    using reference = value_type&;
    // using dim = Dim;

    std::vector<T> v;
    
    using iterator       = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    
    iterator       begin();
    const_iterator begin() const;
    iterator       end();
    const_iterator end()   const;
    
    
    vector()                             : v()        { std::cout << "default\n";}
    vector(std::initializer_list<T> init): v(init)    { std::cout << "init\n"; }
    vector(size_t size)                  : v(size)    { std::cout << "size_t\n";;}
    vector(const vector<T>& other)       : v(other.v) { std::cout << "copy\n";; }
    vector(vector&& other) noexcept      : v(std::move(other.v)) { std::cout << "move\n";; }
    vector(const std::vector<T> other)   : v(other) {std::cout << "std::vector\n";; }
    
    vector(size_t size, const T& value)   : v(size,value) { std::cout << "size_t, value\n";}

    ~vector() { std::cout << "~vector [" << *this << "]"<< std::endl; };

    
    vector<T> operator+(const vector<T>&) const;
    vector<T> operator-(const vector<T>&) const;
    vector<T> operator*(const vector<T>&) const;
    vector<T> operator/(const vector<T>&) const;

    vector<T>& operator+=(const vector<T>&);
    vector<T>& operator-=(const vector<T>&);
    vector<T>& operator*=(const vector<T>&);
    vector<T>& operator/=(const vector<T>&);

    template <typename U>
    vector<T> operator+(const U&) const;
    vector<T> operator-(const T&) const;
    vector<T> operator*(const T&) const;
    vector<T> operator/(const T&) const;

    vector<T>& operator+=(const T&);
    vector<T>& operator-=(const T&);
    vector<T>& operator*=(const T&);
    vector<T>& operator/=(const T&);

    vector<T>& operator=(const vector<T>&) &;
    
    T& operator[](size_t) &;
    const T& operator[](size_t) const &;
    T&& operator[](size_t) &&;
    const T&& operator[](size_t) const &&;

    // auto&& operator[](this auto&& self, size_t pos) {
    //     return std::forward_like<decltype(self)>(self.v[pos])
    // }
    T& at(size_t);
    const T& at(size_t) const;


    
    T dot(const vector<T>&) const;
    size_t size() const;
    vector<size_t> get_size() const;

    void push_back(const T&);
    void push_back(T&&);

    template <typename... Args>
    void emplace_back(Args&&...); // universal reference



};

template <typename T>
size_t vector<T>::size() const {
    return this->v.size();
}

template <typename T>
void vector<T>::push_back(const T& value) {
    std::cout << "lvalue\n"; 
    this->v.push_back(value);
}
template <typename T>
void vector<T>::push_back(T&& value) {
    std::cout << "rvalue\n";
    this->v.push_back(std::move(value));
}


template <typename T>
template <typename... Args>
void vector<T>::emplace_back(Args&&... args) {
    this->v.emplace_back(std::forward<Args>(args)...);
}

template <typename T>
vector<T> vector<T>::operator+(const vector<T>& other) const {
    vector<T> res(this->size());
    typename std::vector<T>::iterator it  = res.begin();
    typename std::vector<T>::iterator end = res.end();
    typename std::vector<T>::const_iterator it1 = this->begin();
    typename std::vector<T>::const_iterator it2 = other.begin();
    while (it != end) {
        *it = (*it1) + (*it2);
        ++it;
        ++it1;
        ++it2;
    }

    return res;
}

template <typename T>
vector<T> vector<T>::operator-(const vector<T>& other) const {
    vector<T> res(this->size());
    typename std::vector<T>::iterator it  = res.begin();
    typename std::vector<T>::iterator end = res.end();
    typename std::vector<T>::const_iterator it1 = this->begin();
    typename std::vector<T>::const_iterator it2 = other.begin();
    while (it != end) {
        *it = (*it1) - (*it2);
        ++it;
        ++it1;
        ++it2;
    }
    return res;
}

template <typename T>
vector<T> vector<T>::operator*(const vector<T>& other) const {
    vector<T> res(this->size());
    typename std::vector<T>::iterator it  = res.begin();
    typename std::vector<T>::iterator end = res.end();
    typename std::vector<T>::const_iterator it1 = this->begin();
    typename std::vector<T>::const_iterator it2 = other.begin();
    while (it != end) {
        *it = (*it1) * (*it2);
        ++it;
        ++it1;
        ++it2;
    }
    return res;
}

template <typename T>
vector<T> vector<T>::operator/(const vector<T>& other) const {
    vector<T> res(this->size());
    typename std::vector<T>::iterator it  = res.begin();
    typename std::vector<T>::iterator end = res.end();
    typename std::vector<T>::const_iterator it1 = this->begin();
    typename std::vector<T>::const_iterator it2 = other.begin();
    while (it != end) {
        *it = (*it1) / (*it2);
        ++it;
        ++it1;
        ++it2;
    }
    return res;
}

template <typename T>
vector<T>& vector<T>::operator+=(const vector<T>& other) {
    typename std::vector<T>::const_iterator it = this->begin();
    typename std::vector<T>::iterator end = this->end();
    typename std::vector<T>::const_iterator it1 = other.begin();
    while (it != end) {
        *it += *it1;
        ++it;
        ++it1;
    }
    return *this;
}

template <typename T>
vector<T>& vector<T>::operator-=(const vector<T>& other) {
    typename std::vector<T>::const_iterator it = this->begin();
    typename std::vector<T>::iterator end = this->end();
    typename std::vector<T>::const_iterator it1 = other.begin();
    while (it != end) {
        *it -= *it1;
        ++it;
        ++it1;
    }
    return *this;
}

template <typename T>
vector<T>& vector<T>::operator*=(const vector<T>& other) {
    typename std::vector<T>::const_iterator it = this->begin();
    typename std::vector<T>::iterator end = this->end();
    typename std::vector<T>::const_iterator it1 = other.begin();
    while (it != end) {
        *it *= *it1;
        ++it;
        ++it1;
    }
    return *this;
}

template <typename T>
vector<T>& vector<T>::operator/=(const vector<T>& other) {
    typename std::vector<T>::const_iterator it = this->begin();
    typename std::vector<T>::iterator end = this->end();
    typename std::vector<T>::const_iterator it1 = other.begin();
    while (it != end) {
        *it /= *it1;
        ++it;
        ++it1;
    }
    return *this;
}

template <typename T>
template <typename U>
vector<T> vector<T>::operator+(const U& other) const {
    vector<T> res(this->size());
    typename std::vector<T>::iterator it  = res.begin();
    typename std::vector<T>::iterator end = res.end();
    typename std::vector<T>::const_iterator it1 = this->begin();
    while (it != end) {
        *it = (*it1) + other;
        ++it;
        ++it1;
    }
    return res;
}

template <typename T>
vector<T> vector<T>::operator-(const T& other) const {
    vector<T> res(this->size());
    typename std::vector<T>::iterator it  = res.begin();
    typename std::vector<T>::iterator end = res.end();
    typename std::vector<T>::const_iterator it1 = this->begin();
    while (it != end) {
        *it = (*it1) - other;
        ++it;
        ++it1;
    }
    return res;
}

template <typename T>
vector<T> vector<T>::operator*(const T& other) const {
    vector<T> res(this->size());
    typename std::vector<T>::iterator it  = res.begin();
    typename std::vector<T>::iterator end = res.end();
    typename std::vector<T>::const_iterator it1 = this->begin();
    while (it != end) {
        *it = (*it1) * other;
        ++it;
        ++it1;
    }
    return res;
}

template <typename T>
vector<T> vector<T>::operator/(const T& other) const {
    vector<T> res(this->size());
    typename std::vector<T>::iterator it  = res.begin();
    typename std::vector<T>::iterator end = res.end();
    typename std::vector<T>::const_iterator it1 = this->begin();
    while (it != end) {
        *it = (*it1) / other;
        ++it;
        ++it1;
    }
    return res;
}

template <typename T>
vector<T>& vector<T>::operator+=(const T& other) {
    typename std::vector<T>::iterator it  = this->begin();
    typename std::vector<T>::iterator end = this->end();
    while (it != end) {
        *it += other;
        ++it;
    }
    return *this;
}

template <typename T>
vector<T>& vector<T>::operator-=(const T& other) {
    typename std::vector<T>::iterator it  = this->begin();
    typename std::vector<T>::iterator end = this->end();
    while (it != end) {
        *it -= other;
        ++it;
    }
    return *this;
}

template <typename T>
vector<T>& vector<T>::operator*=(const T& other) {
    typename std::vector<T>::iterator it  = this->begin();
    typename std::vector<T>::iterator end = this->end();
    while (it != end) {
        *it *= other;
        ++it;
    }
    return *this;
}

template <typename T>
vector<T>& vector<T>::operator/=(const T& other) {
    typename std::vector<T>::iterator it  = this->begin();
    typename std::vector<T>::iterator end = this->end();
    while (it != end) {
        *it /= other;
        ++it;
    }
    return *this;
}

template <typename T>
vector<T> operator+(const T& other, const vector<T>& v) {
    vector<T> res = v;
    v += other;
    return res; // return value optimization
}

template <typename T>
vector<T> operator-(const T& other, const vector<T>& v) {
    vector<T> res(v.size());
    typename std::vector<T>::iterator it  = res.begin();
    typename std::vector<T>::iterator end = res.end();
    typename std::vector<T>::const_iterator it1 = v.begin();
    while (it != end) {
        *it = other - (*it1);
        ++it;
        ++it1;
    }
    return res;
}

template <typename T>
vector<T> operator*(const T& other, const vector<T>& v) {
    vector<T> res = v;
    v *= other;
    return res; // return value optimization
}

template <typename T>
T vector<T>::dot(const vector<T>& other) const {
    T res(0);
    typename std::vector<T>::const_iterator it1 = this->begin();
    typename std::vector<T>::const_iterator end1 = this->end();
    typename std::vector<T>::const_iterator it2 = other.begin();
    while (it1 != end1) {
        res += (*it1) * (*it2);
        ++it1;
        ++it2;
    }
    return res;
}

template <typename T>
T dot(const vector<T>& v1, const vector<T>& v2) {
    T res(0);
    typename std::vector<T>::const_iterator it1 =  v1.begin();
    typename std::vector<T>::const_iterator end1 = v1.end();
    typename std::vector<T>::const_iterator it2 =  v2.begin();
    while (it1 != end1) {
        res += (*it1) * (*it2);
        ++it1;
        ++it2;
    }
    return res;
}

template <typename T>
vector<T>& vector<T>::operator=(const vector<T>& other) & {
    this->v = other.v;
    return *this;
}

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

template<typename T>
std::ostream& operator<<(std::ostream& stream, const vector<T>& vector) {
    for (const T& x: vector) {
        stream << x << " ";
    }
    return stream;
}

template<typename T>
std::ostream& operator<<(std::ostream& stream, const vector<vector<T>>& vector) {
    for (const math::vector<T>& x: vector) {
        stream << x << "\n";
    }
    return stream;
}


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


template <typename T, size_t Dim>
struct pointer_to_vectors {
    using type = std::variant<vector_t<T,Dim>*,
                              typename pointer_to_vectors<T,Dim-1>::type>;
};

template <typename T>
struct pointer_to_vectors<T,1> {
    using type = vector_t<T,1>*;
};

template <typename T, size_t Dim>
using pointer_to_vectors_t = typename pointer_to_vectors<T,Dim>::type;

template <typename T>
struct Debug {
    Debug() = delete;
};


template <typename T>
vector_t<T,1> repmat1(const T& value, size_t n1) {
    return vector_t<T,1>(n1, value);
}

template <typename T>
vector_t<T,2> repmat2(const T& value, size_t n1, size_t n2) {
    return vector_t<T,2>(n1, vector_t<T,1>(n2, value));
}

template <typename T>
vector_t<T,3> repmat3(const T& value, size_t n1, size_t n2, size_t n3) {
    return vector_t<T,3>(n1, vector_t<T,2>(n2, vector_t<T,1>(n3, value)));
}

template <typename T>
vector_t<T,1> zeros1(size_t n1) {
    return vector_t<T,1>(n1);
}

template <typename T>
vector_t<T,2> zeros2(size_t n1, size_t n2) {
    return vector_t<T,2>(n2, zeros1<T>(n1));
}

template <typename T>
vector_t<T,2> zeros3(size_t n1, size_t n2, size_t n3) {
    return vector_t<T,2>(n2, zeros2<T>(n2,n3));
}

template <typename T, size_t Dim>
vector<size_t> size(const vector_t<T,Dim>& vec) {
    vector<size_t> sz(Dim);
    pointer_to_vectors_t<T,Dim> any_vec = &vec;
    
    sz[0] = any_vec->size();
    for (size_t i = 1; i < Dim; ++i) {
        any_vec = &(*any_vec)[0];
        sz[i] = std::get<vector_t<T,Dim>>(any_vec)->size();
    }
    return sz;
}



} // namespace math


/*
int main() {
    using namespace math;
    vector_t<int,1> v = repmat1(0,5);
    vector_t<int,2> vv = repmat2(1,5,3);
    // vector_t<int,3> vvv = repmat3(2,2,3,4);
    // vector_t<double,2> v2 = zeros2<double>(3,4);
    // std::cout << size<int,1>(v);
    pointer_to_vectors_t<int,3> ptr;
    pointer_to_vectors_t<int,3>::ty
    
    using math::pointer_to_vectors_t<int, 3Ui64> = std::variant<math::vector_t<int, 3Ui64> *,
                                                   std::variant<math::vector_t<int, 2Ui64> *,
                                                                math::vector_t<int, 1Ui64> *>> 

    ptr = &v;
    std::cout << (std::get<vector<int>*>(ptr))->size();
    
    
    // Debug<decltype(ptr)>();


    // Debug<decltype(vv)>();
    // multidim_vector<int> v;
    // vector<vector<vector<int>>> vv(2, vector<vector<int>>(3, vector<int>(4)));
    // std::cout << vv;
    
    // vv.emplace_back(5);
    // vv.push_back(vector<int>(5));
    // std::cout << (vv+1) << "\n";
    // for (vector<int>& v : vv) {
    //     std::cout << v + 1<< '\n';
    // }

    // vector<double> v1{1,2,3};
    // vector<double> v2{4,5,6};
    // vector<double> v3 = vector<double>(vector<double>(v1+v2));
    // v3.push_back(v1[0]);
    // v3.push_back(33);
    // std::cout << v3 << '\n' << v1 << '\n' << v2 << '\n';
    
    
    // vector<double> v4 = (v3 + 5.0);
    
    // v1 = v2;
    // vector<double> v3 = v1/2 * v2;
    // int x = v1.dot(v2);
    // std::cout << x << " ";
    // v1.push_back(10.);
    // std::cout << v1 + v2;
    std::cout << "\n" << "end main\n";
}*/