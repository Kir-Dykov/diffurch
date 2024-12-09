#pragma once

#include <math.h>
#include <limits>
#include <array>
#include <iostream>

using namespace std;


template<size_t N>
using Vec = typename conditional<N == 1, double, array<double, N>>::type;

template<size_t N>
using VecArg = typename conditional<N == 1, Vec<1>, const Vec<N>&>::type;

template <size_t N, size_t M>
using VecMap = function<Vec<M>(VecArg<N>)>;

template <size_t N1, size_t N2, size_t M>
using VecMap2 = function<Vec<M>(VecArg<N1>, VecArg<N2>)>;



/**
 * @brief Element-wise sum of two arrays of the same size.
 */
template <typename T, std::size_t N>
std::array<T, N> operator+(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
}

/**
 * @brief Element-wise difference of two arrays of the same size.
 */
template <typename T, std::size_t N>
std::array<T, N> operator-(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
}

/**
 * @brief Product of an array and scalar.
 */
template <typename T, std::size_t N>
std::array<T, N> operator*(const std::array<T, N>& lhs, T rhs) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] * rhs;
    }
    return result;
}

/**
 * @brief Product of a scalar and array.
 */
template <typename T, std::size_t N>
std::array<T, N> operator*(T lhs, const std::array<T, N>& rhs) {
    std::array<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = lhs * rhs[i];
    }
    return result;
}


/**
 * @brief Divide an array by a scalar.
 */
template <typename T, std::size_t N>
std::array<T, N> operator/(const std::array<T, N>& lhs, T rhs) {
    return lhs * (1/rhs); // Is it more efficient to multiply by reciprocal that is calculated once, than perform a multiple divisions?
}



/**
 * @brief Dot product of two arrays.
 */
template <typename T, std::size_t N>
T operator*(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    T result;
    for (std::size_t i = 0; i < N; ++i) {
        result += lhs[i] * rhs[i];
    }
    return result;
}

template <typename ContainerL, typename ContainerR>
decltype(ContainerL{}[0] * ContainerR{}[0]) dot(const ContainerL& lhs, const ContainerR& rhs, size_t size) {
    decltype(lhs[0] * rhs[0]) result{};
    for (std::size_t i = 0; i < size; ++i) {
        result = result + lhs[i] * rhs[i];
    }
    return result;
}



template <typename T, std::size_t N, std::size_t M>
std::array<T, N + M> concatenate(const std::array<T, N>& arr1, const std::array<T, M>& arr2) {
    std::array<T, N + M> result;
    std::copy(arr1.begin(), arr1.end(), result.begin());
    std::copy(arr2.begin(), arr2.end(), result.begin() + N);
    return result;
}

inline void VecCopy(double from, double& to, size_t i = 0) {
    to = from;
}

template <size_t n, typename T>
inline void VecCopy(T from, array<T, n>& to, size_t i = 0) {
    to[i] = from;
}

template <size_t n, size_t m, typename T>
inline void VecCopy(const array<T, n>& from, array<T, m>& to, size_t i = 0) {
    copy(from.begin(), from.end(), to.begin() + i);
}

template <size_t n, typename T>
inline void VecCopy(const array<T, n>& from, vector<T>& to, size_t i = 0) {
    copy(from.begin(), from.end(), to.begin() + i);
}

inline void VecCopy(double from, vector<double>& to, size_t i = 0) {
    to[i] = from;
}


template <size_t n, typename T>
inline double norm(const array<T, n>& v) {
    return sqrt(v*v);
}

inline double norm(double v) {
    return abs(v);
}