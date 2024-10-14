#pragma once

#include "vector.hpp"
#include <cmath>

namespace math {

template <typename It1, typename It2>
requires std::same_as<typename It1::value_type
                    , typename It2::value_type>
class Slice {
public:
    
    using value_type = typename It1::value_type;
    using reference = value_type&;
    static const bool is_vector_like = true;
    using basic_value_type = basic_value_type_t<value_type>;
    using it_difference_type = std::common_type_t<typename It1::difference_type,
                                                  typename It2::difference_type>;
    
private:

    template <bool IsConst, typename It>
    class base_iterator;

public:
    using iterator = base_iterator<false,
                            std::common_type_t<It1,It2>>;
    using const_iterator = base_iterator<true,
                            std::common_type_t<It1,It2>>;

private:
    It1 from;
    It2 to;
    it_difference_type step_;
    size_t sz;




public:
    
    Slice() = default;
    Slice(const It1& from, const It2& to, it_difference_type step = 1);
    void new_slice(const It1& from, const It2& to, it_difference_type step = 1);
    
    size_t size() const;
    it_difference_type step() const;

    constexpr iterator       begin();
    constexpr const_iterator begin() const;
    constexpr const_iterator cbegin() const;
    constexpr iterator       end();
    constexpr const_iterator end() const;
    constexpr const_iterator cend() const; 

    template <ArithmeticVectorsLike<vector<value_type>> Vec>
    Slice<It1,It2>& operator=(const Vec& other) &;
    
    reference operator[](size_t pos);
    const reference operator[](size_t pos) const;

    
};



template <typename It1, typename It2>
Slice<It1,It2>::Slice(const It1& from, const It2& to, it_difference_type step)
        : from(from), to(to), step_(step), sz(std::ceil(double(to-from)/step)) {};

template <typename It1, typename It2>
void Slice<It1,It2>::new_slice(const It1& from_new, const It2& to_new, it_difference_type step_new) {
    from = from_new;
    to = to_new;
    step_ = step_new;
    sz = std::ceil(double(to-from)/step_);
}

template <typename It1, typename It2>
size_t Slice<It1,It2>::size() const {
    return sz;
}

template <typename It1, typename It2>
typename Slice<It1,It2>::it_difference_type  Slice<It1,It2>::step() const {
    return step_;
}

template <typename It1, typename It2>
constexpr typename Slice<It1, It2>::iterator Slice<It1, It2>::begin() {
    return {from,step_};
}

template <typename It1, typename It2>
constexpr typename Slice<It1, It2>::const_iterator Slice<It1, It2>::begin() const {
    return {from,step_};
}

template <typename It1, typename It2>
constexpr typename Slice<It1, It2>::const_iterator Slice<It1, It2>::cbegin() const {
    return {from,step_};
}

template <typename It1, typename It2>
constexpr typename Slice<It1, It2>::iterator Slice<It1, It2>::end() {
    return {to,step_};
}

template <typename It1, typename It2>
constexpr typename Slice<It1, It2>::const_iterator Slice<It1, It2>::end() const {
    return {to,step_};
}

template <typename It1, typename It2>
constexpr typename Slice<It1, It2>::const_iterator Slice<It1, It2>::cend() const {
    return {to,step_};
}

template <typename It1, typename It2>
template <ArithmeticVectorsLike<vector<typename Slice<It1,It2>::value_type>> Vec>
Slice<It1,It2>& Slice<It1, It2>::operator=(const Vec &other) & {
    auto it = Slice<It1,It2>::begin(), end = Slice<It1,It2>::end();
    auto it1 = other.begin();
    while (it < end) {
        *it = *it1;
        ++it;   ++it1;
    }
    return *this;
}

template <typename It1, typename It2>
auto Slice<It1, It2>::operator[](size_t pos)
        -> reference {
    return from[pos*step_];
}


template <typename It1, typename It2>
auto Slice<It1, It2>::operator[](size_t pos) const
        -> const reference {
    return from[pos*step_];
}

template <typename It1, typename It2>
template <bool IsConst, typename It>
class Slice<It1,It2>::base_iterator{
    friend class Slice<It1,It2>;
public:
    using internal_value_type = std::remove_pointer_t<typename It::pointer>;
    using pointer = std::conditional_t<IsConst, 
                                    const internal_value_type*,
                                          internal_value_type*>;
    using reference = std::conditional_t<IsConst, 
                                    const internal_value_type&,
                                          internal_value_type&>;
    using value_type = typename It::value_type;
    using difference_type = typename It::difference_type;
    using iterator_category = std::bidirectional_iterator_tag;
    
private:
    It it;
    difference_type step_;
    
    base_iterator(const It& it, difference_type step_): it(it), step_(step_) {};

public:
    base_iterator() = default;
    base_iterator(const base_iterator&) = default;
    
    operator base_iterator<true,It>() const {
        return {it,step_};
    };

    base_iterator& operator=(const base_iterator&) = default;

    reference operator*() const {
        return *it;
    };

    pointer operator->() const {
        return operator->(it);
    };

    reference operator[](difference_type n) const {
        return it[n*step_];
    }

    base_iterator& operator++() {
        it += step_;
        return *this;
    };
    base_iterator& operator++(int) {
        base_iterator copy = *this;
        it += step_;
        return copy;
    };
    base_iterator& operator+=(difference_type n) {
        it += step_*n;
        return *this;
    };
    base_iterator operator+(difference_type n) const {
        base_iterator res = *this;
        res += step_*n;
        return res;

    };
    base_iterator& operator-=(difference_type n) {
        it -= step_*n;
        return *this;
    };
    base_iterator operator-(difference_type n) const {
        base_iterator res = *this;
        res -= step_*n;
        return res;
    };

    difference_type operator-(const base_iterator& other) const {
        return it - other.it;
    }

    difference_type operator-(const It& other) const {
        return it - other;
    }

    bool operator<(const base_iterator& other) const {
        return it < other.it;
    };
    bool operator<(const It& other) const {
        return it < other;
    };
    bool operator>(const base_iterator& other) const {
        return other.it < it;
    };
    bool operator>(const It& other) const {
        return other < it;
    };
    
    bool operator<=(const base_iterator& other) const {
        return it <= other.it;
    }
    bool operator<=(const It& other) const {
        return it <= other;
    }

    bool operator==(const base_iterator& other) const = delete;
    bool operator==(const It& other) const = delete;
};



} // namespace math