#pragma once
#include <tuple>

using namespace std;

template <typename T>
struct tuple_add_const_s {};

template <typename... Ts>
struct tuple_add_const_s<tuple<Ts...>> {
    using type = std::tuple<typename std::add_lvalue_reference<typename std::add_const<typename std::remove_reference<Ts>::type>::type>::type...>;
};

template <typename T>
struct tuple_decay_s {};

template <typename... Ts>
struct tuple_decay_s<tuple<Ts...>> {
    using type = std::tuple<typename std::decay<Ts>::type...>;
};

template <typename T>
using tuple_add_const = typename tuple_add_const_s<T>::type;

template <typename T>
using tuple_decay = typename tuple_decay_s<T>::type;

template <typename... Ts1, typename... Ts2, size_t... I>
void copy_tuple_implementation(const tuple<Ts1...>& src, tuple<Ts2...>& dest, std::index_sequence<I...>) {
    ((std::get<I>(dest) = std::get<I>(src)), ...);
}

template <typename... Ts1, typename... Ts2>
void copy_tuple(const tuple<Ts1...>& src, tuple<Ts2...>& dest) {
    copy_tuple_implementation(src, dest, std::make_index_sequence<sizeof...(Ts1)>{});
};


template <typename> struct is_tuple: std::false_type {};
template <typename ...T> struct is_tuple<std::tuple<T...>>: std::true_type {};