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

/**
 * @brief like `linspace`, but the values are spaced equally on a logarithmic scale, so the points are more dence near zero.
  * 
 * Generates a vector<T> of `size`  values that are equally spaced on a log-scale in range from `lo` to `hi`. Works as `np.exp(np.linspace(np.log(lo), np.log(hi), size))` in the python package numpy.
 * @param[in] lo  lower bound. Must be positive.
 * @param[in] hi  higher bound. Must be positive.
 * @param[in] size the size of the array.
 * @return vector<T> of size `size`.
 */

template <typename T>
vector<T> expspace(T lo, T hi, size_t size) {
	std::vector<T> v(size);
	for (size_t i = 0; i < size; i++) {
		v[i] = exp( log(lo) + (T(i) / T(size - 1)) * (log(hi) - log(lo)) );
	}
	return v;
}

// template <typename TT> TT sign(TT val) {
// 	return (TT(0) < val) - (val < TT(0));
// }