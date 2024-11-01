#pragma once

using float32 = float;
using float64 = double;
using float128 = long double;

// options: float, double, long double, etc
#ifndef REAL
#define REAL double
#endif

using Real = REAL;


constexpr Real operator"" r(long double value) {
    return static_cast<Real>(value);
}
constexpr Real operator"" r(unsigned long long value) {
    return static_cast<Real>(value);
}