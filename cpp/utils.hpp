#pragma once

#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>

using namespace std;


vector<double> linspace(double lo, double hi, size_t size) {
	std::vector<double> v(size);
	for (size_t i = 0; i < size; i++) {
		v[i] = lo + (double(i) / double(size - 1)) * (hi - lo);
	}
	return v;
}

vector<double> expspace(double lo, double hi, size_t size) {
	std::vector<double> v(size);
	for (size_t i = 0; i < size; i++) {
		v[i] = exp( log(lo) + (double(i) / double(size - 1)) * (log(hi) - log(lo)) );
	}
	return v;
}

template <typename TT> TT sign(TT val) {
	return (TT(0) < val) - (val < TT(0));
}