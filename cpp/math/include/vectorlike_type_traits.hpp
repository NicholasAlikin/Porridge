/*Vector like type traits: concepts, meta-functions*/
#include <concepts>
#include <type_traits>


namespace math {




template <typename T>
concept VectorLike = requires() {
    std::remove_reference_t<T>::is_vector_like;
};

template <typename T>
concept NotVectorLike = (!VectorLike<T>);

/*Iterator of VectorLike*/
template <typename T>
concept VectorLikeIterator = requires(T it, T other) {
    ++it;
    *it;
    it < other;
};
/*Iterator of VectorLike on VectorLike*/
template <typename T>
concept VectorLikeIteratorOn = VectorLikeIterator<T>
&& requires(T it) {
    {*it} -> VectorLike;
};

template <typename T>
concept NotVectorLikeIterator = (!VectorLikeIterator<T>);

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

namespace detail {

template <typename T>
struct basic_value_type_helper {
    using type = T;
};

template <VectorLike vec>
struct basic_value_type_helper<vec> {
    using type = typename basic_value_type_helper<typename vec::value_type>::type;
};

} // namespace detail

template <typename T>
using basic_value_type_t = typename detail::basic_value_type_helper<
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




template <typename ...Iterators>
struct is_arithmetic_iterators {
    static const bool value = false;
};


template <typename ...Iterators>
requires requires(Iterators... types) {
    (*types + ...);
}
struct is_arithmetic_iterators<Iterators...> {
    static const bool value = true;
};
template <typename ...Iterators>
static const bool is_arithmetic_iterators_v = is_arithmetic_iterators<Iterators...>::value;

template <typename ...Iterators>
concept ArithmeticIterators = is_arithmetic_iterators_v<Iterators...>; 

} // namespace math 