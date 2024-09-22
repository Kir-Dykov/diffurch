#pragma once

#include <math.h>
#include <limits>
#include <array>
#include <iostream>

using namespace std;


// template<size_t N>
// using Vec = array<double, N>;

template<size_t N>
using Vec = conditional<N == 1, double, array<double, N>>::type;


// template <size_t N>
// using VecScalarFunc = function<double(const Vec<N>&)>;

// template <size_t N>
// using VecFunc = function<Vec<N>(double)>;

// template <size_t N>
// using VecFunc = function<Vec<N>(double)>;

template <size_t N, size_t M>
using VecMap = function<Vec<M>(const Vec<N>&)>;

template <size_t N1, size_t N2, size_t M>
using VecMap2 = function<Vec<M>(const Vec<N1>&, const Vec<N2>&)>;


// I Want to have Vec<1> === double
// the only inconvinience I see with it
// is that if I want to compose a vector argument
// consisting of many doubles i need to have a copy function with offset:


void VecCopy(double from, double& to, size_t i = 0) {
    to = from;
}

template <size_t n>
void VecCopy(double from, array<double, n>& to, size_t i = 0) {
    to[i] = from;
}

template <size_t n, size_t m>
void VecCopy(array<double, n> from, array<double, m> to, size_t i = 0) {
    copy(from.begin(), from.end(), to.begin() + i);
}


template <size_t n>
void VecCopy(array<double, n> from, vector<double> to, size_t i = 0) {
    copy(from.begin(), from.end(), to.begin() + i);
}

void VecCopy(double from, vector<double>& to, size_t i = 0) {
    to[i] = from;
}

// template <size_t n, size_t k = 0>
// class VecScalarFuncC {
//     array<VecScalarFunc<n>, k+1> fs;
    
// public:
//     VecScalarFuncC(const array<VecScalarFunc<n>, k+1>& fs) : fs(fs) {};
    
//     inline double operator()(Vec<n> t) {
//         return fs[0](t);
//     }
    
//     const VecFunc<n>& operator[](size_t i) {
//         return fs[i];
//     };
// };

// template <size_t n>
// class VecFuncC<n, 0> {
//     VecFunc<n> f;
    
// public:
//     VecFuncC(const VecFunc<n>& f) : f(f) {};
    
//     inline Vec<n> operator()(double t) {
//         return f(t);
//     }
    
//     // added for completeness only, here instead of func[0](t) you should use func(t)
//     const VecFunc<n>& operator[](size_t i) { 
//         if (i > 0) throw("Index is out of range in VecFuncC.");
//         return f;
//     };
// };





// template <size_t n, size_t k = 0>
// class VecFuncC {
//     array<VecFunc<n>, k+1> fs;
    
// public:
//     VecFuncC(const array<VecFunc<n>, k+1>& fs) : fs(fs) {};
    
//     inline Vec<n> operator()(double t) {
//         return fs[0](t);
//     }
    
//     const VecFunc<n>& operator[](size_t i) {
//         return fs[i];
//     };
// };

// template <size_t n>
// class VecFuncC<n, 0> {
//     VecFunc<n> f;
    
// public:
//     VecFuncC(const VecFunc<n>& f) : f(f) {};
    
//     inline Vec<n> operator()(double t) {
//         return f(t);
//     }
    
//     // added for completeness only, here instead of func[0](t) you should use func(t)
//     const VecFunc<n>& operator[](size_t i) { 
//         if (i > 0) throw("Index is out of range in VecFuncC.");
//         return f;
//     };
// };









/**
 * @class MyClass
 * @brief This is a simple class example.
 * 
 * MyClass is used to demonstrate how Doxygen processes class
 * documentation in C++.
 */
template <size_t n, size_t m, size_t k = 0>
class VecMapC {
    VecMap<n, m> f;
    array<VecMap2<n,n,m>, k> df;
    
public:
    
    VecMapC () {};
    
    template <size_t ZERO = 0, typename = std::enable_if_t<k == 0 && ZERO == 0>>
    VecMapC(const VecMap<n, m>& f) : f(f) {};
    
    VecMapC(const VecMap<n, m>& f, const array<VecMap2<n,n,m>, k>& df) : f(f), df(df) {};
    
    inline Vec<m> operator()(const Vec<n>& x) const {
        return f(x);
    }
    
    template <size_t derivative = 0, typename = std::enable_if_t<derivative == 0>>
    inline Vec<m> eval(const Vec<n>& x) const {
        return f(x);
    }
    
    template <size_t derivative, typename = std::enable_if_t<(derivative > 0 && derivative < k)>>
    inline Vec<m> eval(const Vec<n>& x, const Vec<n>& delta_x) const {
        return df[derivative-1](x, delta_x);
    }
};


template <size_t m, size_t k = 0>
class VecFuncC {
    VecMap<1, m> f;
    array<VecMap<1, m>, k> df;
    
public:
    
    VecFuncC () {};
    
    template <size_t ZERO = 0, typename = std::enable_if_t<k == 0 && ZERO == 0>>
    VecFuncC(const VecMap<1, m>& f) : f(f) {};
    
    VecFuncC(const VecMap<1, m>& f, const array<VecMap<1, m>, k>& df) : f(f), df(df) {};
    
    inline Vec<m> operator()(double x) {
        return f(x);
    }
    
    template <size_t derivative = 0>
    inline std::enable_if_t<derivative == 0, Vec<m>> eval(double x) {
        return f(x);
    }
    
    template <size_t derivative = 0>
    inline std::enable_if_t<(derivative > 0 && derivative < k), Vec<m>> eval(double x) {
        return df[derivative-1](x);
    }
    
};






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
    decltype(lhs[0] * rhs[0]) result{};
    for (std::size_t i = 0; i < size; ++i) {
        result = result + lhs[i] * rhs[i];
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









// template<bool is_constant>
// class Delay;

// template<>
// class Delay<true> {
//     double tau;
//     Delay(double tau) : tau(tau) {};
//     inline double operator()(double t) const {
//         return tau;
//     }
// }

// template<>
// class Delay<false> {
//     function<double(double)> tau;
//     Delay(function<double(double)> tau) : tau(tau) {};
//     inline double operator()(double t) const {
//         return tau(t);
//     }
// }
// // HOW WOULD DiscoQue discontinuity propagation work???