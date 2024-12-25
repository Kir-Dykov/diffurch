#pragma once
#include <tuple>
#include <vector>
#include <utility>

using namespace std;

template <typename... Types, size_t... Indices>
void push_back_implementation(tuple<vector<Types>...>& vectors, const tuple<Types...>& values, std::index_sequence<Indices...>) {
     ((std::get<Indices>(vectors).push_back(std::get<Indices>(values))), ...);
}
template <typename... Types>
void push_back(tuple<vector<Types>...>& vectors, 
               const tuple<Types...>& values) {
    push_back_implementation(vectors, values, std::index_sequence_for<Types...>{});
}


template <typename Func, typename... Types>
void for_each(std::tuple<Types...>& tpl, Func&& func) {
    std::apply(
        [&func](auto&... elements) {
            (func(elements), ...); // Fold expression to call func for each element
        },
        tpl
    );
}


// Helper to apply a function to each element of a tuple and create a new tuple
template <typename Func, typename Tuple, std::size_t... Indices>
auto transform_tuple_impl(const Tuple& tpl, Func&& func, std::index_sequence<Indices...>) {
    return std::make_tuple(func(std::get<Indices>(tpl))...); // Apply function and pack results
}

// General function to transform a tuple
template <typename Func, typename... Types>
auto transform(const std::tuple<Types...>& tpl, Func&& func) {
    return transform_tuple_impl(
        tpl, 
        std::forward<Func>(func), 
        std::make_index_sequence<sizeof...(Types)>{}
    );
}



template <typename Condition, typename Tuple>
struct filter_types {}

template <typename Condition>
struct filter_types<Condition, tuple<>> {
    using type = std::tuple<>;  // Return empty tuple
}


template <typename Condition, typename Head, typename... Tail>
struct filter_types<Condition, tuple<Head, Tail...>> {
private:
    using tail_type = typename filter_types<Condition, Tail...>::type;

public:
    using type = std::conditional_t<
        Condition::template value<Head>,  // If condition is satisfied
            decltype(std::tuple_cat(std::tuple<Head>{}, tail_type)),  // Include Head
            tail_type  // Skip Head
    >;
};

// Helper alias for cleaner syntax
template <typename Condition, typename... Ts>
using filter_types = typename filter_types<Condition, Ts...>::type;

template <typename Base>
struct inherits {
    template <typename T>
    static constexpr bool value = std::is_base_of_v<Base, T>;
};


// template <typename T>
// struct tuple_add_const_s {};

// template <typename... Ts>
// struct tuple_add_const_s<tuple<Ts...>> {
//     using type = std::tuple<typename std::add_lvalue_reference<typename std::add_const<typename std::remove_reference<Ts>::type>::type>::type...>;
// };

// template <typename T>
// struct tuple_decay_s {};

// template <typename... Ts>
// struct tuple_decay_s<tuple<Ts...>> {
//     using type = std::tuple<typename std::decay<Ts>::type...>;
// };

// template <typename T>
// using tuple_add_const = typename tuple_add_const_s<T>::type;

// template <typename T>
// using tuple_decay = typename tuple_decay_s<T>::type;

// template <typename... Ts1, typename... Ts2, size_t... I>
// void copy_tuple_implementation(const tuple<Ts1...>& src, tuple<Ts2...>& dest, std::index_sequence<I...>) {
//     ((std::get<I>(dest) = std::get<I>(src)), ...);
// }

// template <typename... Ts1, typename... Ts2>
// void copy_tuple(const tuple<Ts1...>& src, tuple<Ts2...>& dest) {
//     copy_tuple_implementation(src, dest, std::make_index_sequence<sizeof...(Ts1)>{});
// };


// template <typename> struct is_tuple: std::false_type {};
// template <typename ...T> struct is_tuple<std::tuple<T...>>: std::true_type {};


