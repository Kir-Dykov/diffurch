#pragma once

#include <math.h>
#include <limits>
#include <array>
#include <iostream>

using namespace std;


template<size_t N>
using Vec = array<double, N>;

// template<size_t N>
// using Vec = conditional<N == 1, double, array<double, N>>::type;


template <size_t N>
using VecScalarFunc = function<double(const Vec<N>&)>;

template <size_t N>
using VecFunc = function<Vec<N>(double)>;


// template <size_t N>
// using VecFunc = function<Vec<N>(double)>;

template <size_t N, size_t M>
using VecMap = function<Vec<M>(const Vec<N>&)>;

template <size_t N1, size_t N2, size_t M>
using VecMap2 = function<Vec<M>(const Vec<N1>&, const Vec<N2>&)>;


// template <size_t N, size_t M>
// using ArrayFunc = function<array<double, M>(const array<double, N>&)>;

// template <size_t N>
// using ArrayFuncScalar = function<double(const array<double, N>&)>;

// template <size_t N1, size_t N2, size_t M>
// using ArrayFunc2 = function<array<double, M>(const array<double, N1>&,const array<double, N2>&)>;

// print Array
template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {
    os << '[';
    for (std::size_t i = 0; i < N; ++i) {
        os << arr[i];
        if (i < N - 1) {
            os << ", ";
        }
    }
    os << ']';
    return os;
}


// Array + Array
template <typename T, std::size_t N>
std::array<T, N> operator+(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
}

// Array - Array
template <typename T, std::size_t N>
std::array<T, N> operator-(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
}

// Array * Scalar
template <typename T, std::size_t N>
std::array<T, N> operator*(const std::array<T, N>& lhs, T rhs) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] * rhs;
    }
    return result;
}



// Array / Scalar
template <typename T, std::size_t N>
std::array<T, N> operator/(const std::array<T, N>& lhs, T rhs) {
    return lhs * (1/rhs);
}

// Array * Scalar
template <typename T, std::size_t N>
std::array<T, N> operator*(T lhs, const std::array<T, N>& rhs) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs * rhs[i];
    }
    return result;
}

// dot product :: Array * Array
template <typename T, std::size_t N>
T operator*(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    T result;
    for (std::size_t i = 0; i < N; ++i) {
        result += lhs[i] * rhs[i];
    }
    return result;
}

// dot product :: Vector * Vector
template <typename ContainerL, typename ContainerR>
decltype(ContainerL{}[0] * ContainerR{}[0]) dot(const ContainerL& lhs, const ContainerR& rhs, size_t size) {
    decltype(lhs[0] * rhs[0]) result;
    for (std::size_t i = 0; i < size; ++i) {
        result += lhs[i] * rhs[i];
    }
    return result;
}




template <typename T, std::size_t N, std::size_t M>
std::array<T, N + M> concatenate(const std::array<T, N>& arr1, const std::array<T, M>& arr2) {
    std::array<T, N + M> result;

    // Copy the first array into the result
    std::copy(arr1.begin(), arr1.end(), result.begin());

    // Copy the second array into the result, starting after the first array
    std::copy(arr2.begin(), arr2.end(), result.begin() + N);

    return result;
}