#pragma once

#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>
#include <assert.h>

using namespace std;



/*
uint8 : number of arrays

char : type
    'd':double, 'f':float, 'b'/'c':byte/char, 'i'/'u':int, 'I'/'U':long int
uint32 :           dimension
uint32*dimension : shape
T*size : data

*/


template <typename T>
struct type_char {};
template<>
struct type_char<char> { static const char value = 'c'; };
template<>
struct type_char<unsigned char> { static const char value = 'b'; };
template<>
// struct type_char<short> { static const char value = 's'; };
// template<>
struct type_char<int> { static const char value = 'i'; };
template<>
struct type_char<long long> { static const char value = 'I'; };
template<>
struct type_char<unsigned int> { static const char value = 'u'; };
template<>
struct type_char<unsigned long long> { static const char value = 'U'; };
template<>
struct type_char<float> { static const char value = 'f'; };
template<>
struct type_char<double> { static const char value = 'd'; };
template<>
struct type_char<long double> { static const char value = 'F'; };
template <typename T>
using type_char_v = type_char<T>::value;



template<typename T>           struct ArrayType                   { using type = T; };
template<typename T>           struct ArrayType<std::vector<T>>   { using type = typename ArrayType<T>::type; };
template<typename T, size_t n> struct ArrayType<std::array<T, n>> { using type = typename ArrayType<T>::type; };

template<typename T>           struct ArraySpec                   { static const int dimension = 0;                           using type = T;  };
template<typename T>           struct ArraySpec<std::vector<T>>   { static const int dimension = 1 + ArraySpec<T>::dimension; using type = typename ArraySpec<T>::type; };
template<typename T, size_t n> struct ArraySpec<std::array<T, n>> { static const int dimension = 1 + ArraySpec<T>::dimension; using type = typename ArraySpec<T>::type; };

template<typename T>
vector<size_t> shape_of(T arr) {
    if constexpr (ArraySpec<T>::dimension == 0) {
        return {}; // not an array
    } else {
        vector<size_t> result {arr.size()};
        assert(result[0] > 0 && "Array must not be empty.");
        vector<size_t> append = shape_of(arr[0]);
        result.insert(result.end(), append.begin(), append.end());
        return result;
    }
}

template <typename T, typename U>
void flatten_into(const T& vec, vector<U>& result) {
    if constexpr (ArraySpec<T>::dimension == 1) {
        result.insert(result.end(), vec.begin(), vec.end());
    } else {
        for (const auto& element : vec) {
            flatten_into(element, result);
        }
    }
}


void write_arrays(ofstream& file) {
    return;
}
template <typename Array, typename... Rest>
void write_arrays(ofstream& file, const Array& array, Rest... rest) {

    using array_type = typename ArraySpec<Array>::type;

    vector<size_t> shape = shape_of(array);
    size_t dimension     = shape.size();
    assert(dimension > 0  &&  "Dimension must be at least 1");
    char type_c = type_char<array_type>::value;
    
    file.write((char*)(&type_c),       sizeof(char)              );
    file.write((char*)(&dimension),    sizeof(size_t)            );
    file.write((char*)(shape.data()), sizeof(size_t) * dimension);
    
    vector<array_type> data;
    flatten_into(array, data);
        
    file.write((char*)(data.data()), sizeof(array_type)*data.size());
    
    write_arrays(file, rest...);
}

template <typename... Types>
void save_arrays(string filename, tuple<Types...> arrays) {
    save_arrays(filename, std::get<>(arrays)...)
}

template <typename... Types>
void save_arrays(string filename, Types... arrays) {
    ofstream file(filename, ios::binary); // open file to write to
    assert(bool(file));
    
    size_t number_of_arrays = sizeof...(Types);
    file.write((char*)(&number_of_arrays),    sizeof(size_t));
    
    write_arrays(file, arrays...);
}






















// template <typename T> void save(vector<T> image, string filename) {
    
    
//     ofstream file(filename, ios::binary); // open file to write to
    
//     // cout << "size of T is " << sizeof(T) << endl;
//     if (!file) {
//         cout << "Failed to open the file: " << filename << endl;
//         return;
//     }
    
//     file.write((char*)(&(typeid(T).name()[0])), sizeof(char)); // type of array
    
//     int shape_size = 1;
//     int       size1= image.size();
    
// 	file.write((char*)(&shape_size), sizeof(int)); // dimentionality
    
// 	file.write((char*)(&size1), sizeof(int)); // size
    
//     file.write((char*)(&image[0]), sizeof(T) * size1);    // array itself
// }

// template <typename T> void save(vector<vector<T>> image, string filename) {
    
//     ofstream file(filename, ios::binary); // open file to write to
//      if (!file) {
//         cout << "Failed to open the file: " << filename << endl;
//         return;
//     }
    
//     file.write((char*)(&(typeid(T).name()[0])), sizeof(char)); // type of array
    
    
    
//     int shape_size  = 2;
//     int       size1 = image.size();
//     int       size2 = image[0].size();
    
// 	file.write((char*)(&shape_size), sizeof(int)); // dimentionality
    
// 	file.write((char*)(&size1), sizeof(int)); // size
// 	file.write((char*)(&size2), sizeof(int)); // size
    
//     for (int i = 0; i < size1; i++) {
//         file.write((char*)(&(image[i][0])), sizeof(T) * size2);    // array itself
//     }
// }

// template <typename T> void save(vector<vector<vector<T>>> image, string filename) {
//     ofstream file(filename, ios::binary); // open file to write to
//      if (!file) {
//         cout << "Failed to open the file: " << filename << endl;
//         return;
//     }
    
//     file.write((char*)(&(typeid(T).name()[0])), sizeof(char)); // type of array
    
    
//     int shape_size  = 3;
//     int       size1 = image.size();
//     int       size2 = image[0].size();
//     int       size3 = image[0][0].size();
    
// 	file.write((char*)(&shape_size), sizeof(int)); // dimentionality
    
// 	file.write((char*)(&size1), sizeof(int)); // size
// 	file.write((char*)(&size2), sizeof(int)); // size
// 	file.write((char*)(&size3), sizeof(int)); // size
    
//     for (int i = 0; i < size1; i++) {
//         for (int j = 0; j < size2; j++)
//             file.write((char*)(&(image[i][j][0])), sizeof(T) * size3);    // array itself
//     }
// }
