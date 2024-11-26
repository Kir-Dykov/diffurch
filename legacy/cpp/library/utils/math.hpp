#pragma once

#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>

using namespace std;


/**
 * @brief array of linearly spaced values.
 * 
 * Generates a vector<T> of `size` equally spaced values in range from `lo` to `hi`. Works as a function `linspace` from `python`
 package `numpy`. * 
 * @param[in] lo lower bound.
 * @param[in] hi higher bound.
 * @param[in] size the size of the array.
 * @return vector<T> of size `size`.
 */

template <typename T>
vector<T> linspace(T lo, T hi, size_t size) {
	std::vector<T> v(size);
	for (size_t i = 0; i < size; i++) {
		v[i] = lo + (T(i) / T(size - 1)) * (hi - lo);
	}
	return v;
}


template <typename T>
vector<T> expspace(T lo, T hi, size_t size) {
	std::vector<T> v(size);
	for (size_t i = 0; i < size; i++) {
		v[i] = exp( log(lo) + (T(i) / T(size - 1)) * (log(hi) - log(lo)) );
	}
	return v;
}

template <typename TT> TT sign(TT val) {
	return (TT(0) < val) - (val < TT(0));
}

template <typename TT> int minus_one_to_the(TT val) {
    return 1 - 2 * (val & 1);
}


template <typename R, typename T>
function<R(T)> periodic_continuation(T a, T b, function<R(T)> f) {
    return [a,b,f](T x){
        double p = b - a;
        x = fmod(x - a, p);
        x += (x < 0)*p + a;
        
        return f(x);
    };
}

template <typename T>
T fmod_period(T x, T a, T b) {
    T p = b - a;
    x = fmod(x - a, p);
    return x +  (x < 0)*p + a;
}

template <typename T>
T fmod_positive(T x, T p) {
    x = fmod(x, p);
    return x +  (x < 0)*p;
}