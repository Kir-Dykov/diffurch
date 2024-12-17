#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <tuple>
#include <vector>
#include <type_traits>
#include <iostream>

namespace py = pybind11;

using namespace std;

// Forward declarations
template <typename T>
auto convert_to_numpy(const T& value);
template <typename... Args>
auto convert_to_numpy(const std::tuple<Args...>& tuple);
template <typename T>
auto convert_to_numpy(const std::vector<T>& vec);

// Specialization for std::tuple
template <typename Tuple, std::size_t... Is>
auto convert_to_numpy_tuple_implementation(const Tuple& tuple, std::index_sequence<Is...>) {
    return std::make_tuple(convert_to_numpy(std::get<Is>(tuple))...);
}
template <typename... Args>
auto convert_to_numpy(const std::tuple<Args...>& tuple) {
    return convert_to_numpy_tuple_implementation(tuple, std::index_sequence_for<Args...>{});
}

// Specialization for std::vector
template <typename T>
auto convert_to_numpy(const std::vector<T>& vec) {
    return py::array(vec.size(), vec.data());
}

// Generic conversion (fallback for scalar types or unsupported types)
template <typename T>
auto convert_to_numpy(const T& value) {
    return value;
}
