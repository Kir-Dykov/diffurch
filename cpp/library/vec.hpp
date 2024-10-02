#pragma once

#include <math.h>
#include <limits>
#include <array>
#include <iostream>
#include "real.hpp"

using namespace std;


/**
 * @brief Alias for array<Real, N> if N > 1 and Real if N == 1.
 */
template<size_t N>
using Vec = conditional<N == 1, Real, array<Real, N>>::type;

/**
 * @brief Alias for type of function from Vec<N> to Vec<M>.
   @see Vec
 */
template <size_t N, size_t M>
using VecMap = function<Vec<M>(const Vec<N>&)>;

/**
 * @brief Alias for type of function from Vec<N1> and Vec<N2> to Vec<M>.
 This is useful for defining the high-dimensional derivatives for VecMap's.
   @see Vec
   @see VecMap
 */
template <size_t N1, size_t N2, size_t M>
using VecMap2 = function<Vec<M>(const Vec<N1>&, const Vec<N2>&)>;



/**
 * @brief Class encapsulating the function from `Vec<n>` to `Vec<m>` as well as `k` it's derivatives.
 @tparam n The dimension for arguments of function.
 @tparam m The dimension for values of function.
 @tparam k The number of derivatives provided for that function.
 */
template <size_t n, size_t m, size_t k = 0>
class VecMapC {
    using derivative_type = std::conditional_t<(n > 1), VecMap2<n,n,m>, VecMap<n,m> >;
    VecMap<n, m> f;
    // if the argument of function is scalar, then, instead of differential operator we just have higher order total derivatives
    array<derivative_type, k> df;
    
public:
    
    VecMapC () {};
    
    template <size_t ZERO = 0, typename = std::enable_if_t<k == 0 && ZERO == 0>>
    VecMapC(const VecMap<n, m>& f) : f(f) {};
    
    VecMapC(const VecMap<n, m>& f, const array<derivative_type, k>& df) : f(f), df(df) {};
    
    /**
     * @brief Evaluate function at the value `x`
    */
    inline Vec<m> operator()(const Vec<n>& x) const {
        return f(x);
    }
    
    
    /**
     * @brief Apply function's derivative at the value `x` to a finite displacement 'delta_x'.
     For the case of `derivative_order == 1' this evaluates to `f'(x) * delta_x`, where f'(x) is a jacobian matrix.
     For the case of `derivative_order == 2' this evaluates to `delta_x^T * f''(x) * delta_x`, where `f''(x)` is a matrix of second order partial derivatives.
     @tparam derivative_order The order of differential operator
    */
    
    
    
    
    template <size_t derivative_order = 0>
    inline Vec<m> eval(const Vec<n>& x) const {
        if constexpr (derivative_order == 0) {
            return f(x);
        } else if constexpr (n == 1) {
            return df[derivative_order-1](x);
        } else {
            static_assert(derivative_order > 0 && n > 1, "eval<k> for k > 0 takes two arguments when n > 1");
        }
    }
        
    template <size_t derivative_order>
    inline std::enable_if_t<(derivative_order > 0 && derivative_order < k), Vec<m>> 
    eval(const Vec<n>& x, const Vec<n>& delta_x) const {
        if constexpr (n == 1) {
            return df[derivative_order-1](x) * delta_x;
        } else {
            return df[derivative_order-1](x, delta_x);
        }
    }
    
    // template <size_t derivative_order = 0, typename Container>
    vector<Vec<m>> eval_series(vector<Vec<n>> data) {
        vector<Vec<m>> result(data.size());
        std::transform(data.begin(), data.end(), result.begin(), [this](const Vec<n>& x){return eval(x);});
        return result;
    }
    
    
    
};











/**
 * @brief Operator to print array in the format [1, 2, 3]
 */
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

/**
 * @brief Generalized dot product, usefull to define Runge Kutta method.
 @param[in] size Number of elements of each argument that are used to 
 @see RK_TimeSeries::step_push
 */
template <typename ContainerL, typename ContainerR>
decltype(ContainerL{}[0] * ContainerR{}[0]) dot(const ContainerL& lhs, const ContainerR& rhs, size_t size) {
    decltype(lhs[0] * rhs[0]) result{};
    for (std::size_t i = 0; i < size; ++i) {
        result = result + lhs[i] * rhs[i];
    }
    return result;
}



/**
 * @brief Concatenation of two arrays.
 */
template <typename T, std::size_t N, std::size_t M>
std::array<T, N + M> concatenate(const std::array<T, N>& arr1, const std::array<T, M>& arr2) {
    std::array<T, N + M> result;

    // Copy the first array into the result
    std::copy(arr1.begin(), arr1.end(), result.begin());

    // Copy the second array into the result, starting after the first array
    std::copy(arr2.begin(), arr2.end(), result.begin() + N);

    return result;
}







/**
 * @brief std::copy extension to work with type alias Vec<n>
 @see Vec
 */
void VecCopy(Real from, Real& to, size_t i = 0) {
    to = from;
}

/**
 * @brief std::copy extension to work with type alias Vec<n>
 @see Vec
 */
template <size_t n, typename T>
void VecCopy(T from, array<T, n>& to, size_t i = 0) {
    to[i] = from;
}

/**
 * @brief std::copy extension to work with type alias Vec<n>
 @see Vec
 */
template <size_t n, size_t m, typename T>
void VecCopy(const array<T, n>& from, array<T, m>& to, size_t i = 0) {
    copy(from.begin(), from.end(), to.begin() + i);
}

/**
 * @brief std::copy extension to work with type alias Vec<n>
 @see Vec
 */
template <size_t n, typename T>
void VecCopy(const array<T, n>& from, vector<T>& to, size_t i = 0) {
    copy(from.begin(), from.end(), to.begin() + i);
}

/**
 * @brief std::copy extension to work with type alias Vec<n>
 @see Vec
 */
void VecCopy(Real from, vector<Real>& to, size_t i = 0) {
    to[i] = from;
}



template <size_t n, typename T>
inline Real norm(const array<T, n>& v) {
    return sqrt(v*v);
}

inline Real norm(Real v) {
    return abs(v);
}