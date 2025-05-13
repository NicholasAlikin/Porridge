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
    using const_reference = const value_type&;
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
	Slice(const It1& from, const It2& to, it_difference_type step, size_t sz);
    
    Slice() = default;
    Slice(const It1& from, const It2& to, it_difference_type step = 1);
    void new_slice(const It1& from, const It2& to, it_difference_type step = 1);
    void update_from(const It1& from);
    void update_to(const It2& to);
    void update_from_to(const It1& from, const It2& to);
    
    size_t size() const;
    it_difference_type step() const;

    constexpr iterator       begin();
    constexpr const_iterator begin() const;
    constexpr const_iterator cbegin() const;
    constexpr iterator       end();
    constexpr const_iterator end() const;
    constexpr const_iterator cend() const; 

    // template <ArithmeticVectorsLike<vector<Slice<It1,It2>::value_type>> Vec>
    // requires std::same_as< decltype(Slice<It1,It2>::begin()),
    //                     typename Slice<It1,It2>::iterator >
    // Slice<It1,It2>& operator=(const Vec& other) const &;
    template <ArithmeticVectorsLike<vector<typename Slice<It1,It2>::value_type>> Vec>
    requires requires(Slice<It1,It2> sl) {
        *(sl.begin()) = Slice<It1,It2>::value_type{};
    }
    Slice<It1, It2>& operator=(const Vec& other) & {
        auto it = begin(), end_ = end();
        auto it1 = other.begin();
        while (it < end_) {
            *it = *it1;
            ++it;   ++it1;
        }
        return *this;
    }


    Slice<It1, It2>& operator=(const value_type& other) & {
        auto it = begin(), end_ = end();
        while (it < end_) {
            *it = other;
            ++it;
        }
        return *this;
    }

    template <ArithmeticVectorsLike<vector<typename Slice<It1,It2>::value_type>> Vec>
    requires requires(Slice<It1,It2> sl) {
        *(sl.begin()) = Slice<It1,It2>::value_type{};
    }
    Slice<It1, It2>& operator=(Vec&& other) & {
        auto it = begin(), end_ = end();
        auto it1 = other.begin();
        while (it < end_) {
            *it = std::move(*it1);
            ++it;   ++it1;
        }
        return *this;
    }


    // template <ArithmeticVectorsLike<vector<value_type>> Vec>
    // Slice<It1,It2>& operator=(const Vec& other) const & = delete;
    
    reference operator[](size_t pos);
    const_reference operator[](size_t pos) const;
	
	
	
	operator Slice<
		std::conditional_t<std::is_const_v<std::remove_reference_t<typename It1::reference>>,It1,
			std::conditional_t<
				std::same_as<It1,typename math::vector<typename It1::value_type>::iterator>
				,typename math::vector<typename It1::value_type>::const_iterator
				,typename math::Slice<It1,It1>::const_iterator
			>
		>,
		std::conditional_t<std::is_const_v<std::remove_reference_t<typename It2::reference>>,It2,
			std::conditional_t<
				std::same_as<It2,typename math::vector<typename It2::value_type>::iterator>
				,typename math::vector<typename It2::value_type>::const_iterator
				,typename math::Slice<It2,It2>::const_iterator
			>
		>
	>() const {
		return {from, to, step_, sz};
	}
        // return {static_cast<ConstIt_t<It1>>(from), static_cast<ConstIt_t<It2>>(to), step_, sz};
    
};

template <typename T, size_t Dim=1>
using vector_slice = Slice<typename vector_t<T,Dim>::iterator,
                           typename vector_t<T,Dim>::iterator>;

template <typename T, size_t Dim=1>
using vector_const_slice = Slice<typename vector_t<T,Dim>::const_iterator,
                                 typename vector_t<T,Dim>::const_iterator>;



// template <typename It1, typename It2>
// Slice<It1,It2>::operator Slice<ConstIt_t<It1>,ConstIt_t<It2>>() const {
	// return {from, to, step_, sz};
// }


///////////
// template <typename It>
// concept ConstIterator = std::is_const_v<typename It::reference>;

// template <ConstIterator It>
// struct ConstIt<It> {
	// using type = It;
// };


// // vector iterator
// template <typename It>
// requires std::same_as<It,typename math::vector<typename It::value_type>::iterator>
// struct ConstIt<It> {
	// using type = typename math::vector<typename It::value_type>::const_iterator;
// };
// // Slice iterator
// template <typename It>
// requires std::same_as<It,typename math::Slice<It,It>::iterator>
// struct ConstIt<It> {
	// using type = typename math::Slice<It,It>::const_iterator;
// };


//////////

template <typename It1, typename It2>
Slice<It1,It2>::Slice(const It1& from, const It2& to, it_difference_type step, size_t sz)
        : from(from), to(to), step_(step), sz(sz) {};

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
void Slice<It1,It2>::update_from(const It1& from_new) {
    from = from_new;
}

template <typename It1, typename It2>
void Slice<It1,It2>::update_to(const It2& to_new) {
    to = to_new;
}

template <typename It1, typename It2>
void Slice<It1,It2>::update_from_to(const It1& from_new, const It2& to_new) {
    from = from_new;
    to = to_new;
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

// template <typename It1, typename It2>
// template <ArithmeticVectorsLike<vector<value_type>> Vec>
// inline Slice<It1, It2> &math::Slice<It1, It2>::operator=(const Vec &other) &
// {
//     // TODO: insert return statement here
// }

// template <typename It1, typename It2>
// template <ArithmeticVectorsLike<vector<typename Slice<It1,It2>::value_type>> Vec>
// requires std::same_as< decltype(Slice<It1,It2>::begin()),
//                        typename Slice<It1,It2>::iterator >
// Slice<It1,It2>& Slice<It1,It2>::operator=(const Vec &other) & {
//     auto it = begin(), end_ = end();
//     auto it1 = other.begin();
//     while (it < end_) {
//         *it = *it1;
//         ++it;   ++it1;
//     }
//     return *this;
// }

template <typename It1, typename It2>
auto Slice<It1, It2>::operator[](size_t pos)
        -> reference {
    return from[pos*step_];
}


template <typename It1, typename It2>
auto Slice<It1, It2>::operator[](size_t pos) const
        -> const_reference {
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
        return it.operator->();
    };

    reference operator*() {
        return *it;
    };

    pointer operator->() {
        return it.operator->();
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
    }
    bool operator<(const It& other) const {
        return it < other;
    }
    bool operator>(const base_iterator& other) const {
        return other.it < it;
    }
    bool operator>(const It& other) const {
        return other < it;
    }

    /*operator "!="  is the same as operator "<"
      because step might be != 1
      but std functions with iterators uses only "!=" operator.*/
    bool operator!=(const base_iterator& other) const {
        return it < other.it;
    }
    bool operator!=(const It& other) const {
        return it < other;
    }
    
    bool operator<=(const base_iterator& other) const {
        return it <= other.it;
    }
    bool operator<=(const It& other) const {
        return it <= other;
    }

    bool operator==(const base_iterator& other) const = delete;
    bool operator==(const It& other) const = delete;
};


// template <typename It>
// requires It::
// bool operator<(const typename vector<T>::const_iterator& it1,
//                const typename Slice<typename vector<T>::const_iterator,
//                                     typename vector<T>::const_iterator>::const_iterator& it2)
// {
//     return it2 > it1;
// }


} // namespace math