#pragma once

#include <math.h>
#include <limits>
#include <array>
#include <iostream>
#include "real.hpp"


#include <boost/preprocessor.hpp>

#define VEC_LAMBDA_many(capture, arg_tuple, value_tuple) \
    [BOOST_PP_TUPLE_REM_CTOR(capture)](VecArg<BOOST_PP_TUPLE_SIZE(arg_tuple)> VEC_LAMBDA_many_arguments){ \
        const auto& [BOOST_PP_TUPLE_REM_CTOR(arg_tuple)] = VEC_LAMBDA_many_arguments;\
        return Vec<BOOST_PP_TUPLE_SIZE(value_tuple)>{BOOST_PP_TUPLE_REM_CTOR(value_tuple)}; \
    }\

#define VEC_LAMBDA_one(capture, arg, value_tuple) \
    [BOOST_PP_TUPLE_REM_CTOR(capture)](VecArg<1> arg){ \
        return Vec<BOOST_PP_TUPLE_SIZE(value_tuple)>{BOOST_PP_TUPLE_REM_CTOR(value_tuple)}; \
    }\

#define VEC_LAMBDA(capture, arg_tuple, value_tuple) \
 BOOST_PP_IIF(BOOST_PP_EQUAL(1, BOOST_PP_TUPLE_SIZE(arg_tuple)), \
    (VEC_LAMBDA_one(capture, BOOST_PP_TUPLE_REM_CTOR(arg_tuple), value_tuple)), \
    (VEC_LAMBDA_many(capture,arg_tuple,value_tuple)) )


using namespace std;


template<size_t N>
using Vec = conditional<N == 1, Real, array<Real, N>>::type;

template<size_t N>
using VecArg = conditional<N == 1, Vec<1>, const Vec<N>&>::type;

template <size_t N, size_t M>
using VecMap = function<Vec<M>(VecArg<N>)>;

template <size_t N1, size_t N2, size_t M>
using VecMap2 = function<Vec<M>(VecArg<N1>, VecArg<N2>)>;



template <size_t n, size_t m, size_t k = 0>
class VecMapC {
 public:

    using derivative_type = std::conditional_t<(n > 1), VecMap2<n,n,m>, VecMap<n,m> >;
    
    VecMap<n, m> f;
    
    // if the argument of function is scalar, then, instead of differential operator we just have higher order total derivatives
    array<derivative_type, k> df;
        
    VecMapC () {};
    
    template <typename Callable, size_t ZERO = 0, typename = std::enable_if_t<k == 0 && ZERO == 0>>
    VecMapC(const Callable& f) : f(f) {};
    
    template <typename Callable1, typename Callable2>
    VecMapC(const Callable1& f, const array<Callable2, k>& df) : f(f), df(df) {};
    
    inline Vec<m> operator()(VecArg<n> x) const {
        return f(x);
    }
    
    template <size_t derivative_order = 0>
    inline Vec<m> eval(VecArg<n> x) const {
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
    eval(VecArg<n> x, VecArg<n> delta_x) const {
        if constexpr (n == 1) {
            return df[derivative_order-1](x) * delta_x;
        } else {
            return df[derivative_order-1](x, delta_x);
        }
    }
    
    template <size_t derivative_order = 0>
    vector<Vec<m>> eval_series(const vector<Vec<n>>& data) {
        vector<Vec<m>> result(data.size());
        std::transform(data.begin(), data.end(), result.begin(), [this](const Vec<n>& x){return eval<derivative_order>(x);});
        return result;
    }
    
    
    
};


template <size_t n>
VecMapC<1,n> ConstantVecMapC(Vec<n> x) {
    return [x](Real t){return x;};
}


template <size_t derivatives_n = 0>
VecMapC<1,1,derivatives_n> CosVecMapC(Real A, Real w, Real gamma = 0) {
    array<VecMap<1,1>, derivatives_n> derivatives;
    int k = 1;
    Real W = w;
    while (k <= derivatives_n) {
        derivatives[k-1] = [A, w, gamma, W](Real t){return - A * W * sin(w * t + gamma);};
        k++; W*=w; if (k > derivatives_n) break;
        derivatives[k-1] = [A, w, gamma, W](Real t){return - A * W * cos(w * t + gamma);};
        k++; W*=w; if (k > derivatives_n) break;
        derivatives[k-1] = [A, w, gamma, W](Real t){return + A * W * sin(w * t + gamma);};
        k++; W*=w; if (k > derivatives_n) break;
        derivatives[k-1] = [A, w, gamma, W](Real t){return + A * W * cos(w * t + gamma);};
        k++; W*=w; if (k > derivatives_n) break;
    }
    return VecMapC<1,1,derivatives_n>([A, w, gamma](Real t){return A * cos(w * t + gamma);}, derivatives);
}

VecMapC<1,2> Cos2VecMapC(Real A, Real w) {
    return VecMapC<1,2>([A, w](Real t){return Vec<2>{A * cos(w * t), w * A * sin(w * t)};});
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







inline void VecCopy(Real from, Real& to, size_t i = 0) {
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

inline void VecCopy(Real from, vector<Real>& to, size_t i = 0) {
    to[i] = from;
}





template <size_t n, typename T>
inline Real norm(const array<T, n>& v) {
    return sqrt(v*v);
}

inline Real norm(Real v) {
    return abs(v);
}