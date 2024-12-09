#pragma once

#include <iostream>
#include <vector>
#include <tuple>

// // Helper function to combine elements of Cartesian product
// template <typename T, typename... Rest>
// void cartesian_product_recursion(const std::vector<std::vector<T>>& vectors, 
//                      size_t depth, 
//                      std::vector<T>& current, 
//                      std::vector<std::vector<T>>& result) {
//     if (depth == vectors.size()) {
//         result.push_back(current);
//         return;
//     }
//     for (const auto& val : vectors[depth]) {
//         current.push_back(val);
//         cartesian_product_recursion(vectors, depth + 1, current, result);
//         current.pop_back();
//     }
// }

// // Main Cartesian product function
// template <typename T>
// std::vector<std::vector<T>> cartesian_product(const std::vector<std::vector<T>>& vectors) {
//     std::vector<std::vector<T>> result;
//     if (vectors.empty()) return result;

//     std::vector<T> current;
//     current.reserve(vectors.size());
//     cartesian_product_recursion(vectors, 0, current, result);
//     return result;
// }



// template <typename T>
// std::vector<std::vector<T>> cartesian_product(const std::vector<std::vector<T>>& vectors) {
    
//     size_t n = vectors.size();  
//     vector<size_t> sizes(n);
//     std::transform(vectors.begin(), vectors.end(), sizes.begin(), [](const vector<T>& v){return v.size();});
//     vector<size_t> indexes(n);
    
//     vector<vector<T>> result;
    
//     size_t result_size = 1;
//     for (int i = 0; i < n; i++) {
//         result_size *= sizes[i];
//     }
//     if (result_size == 0) {
//         return result;
//     }
    
//     result.resize(result_size, vector<T>(n));
//     int j = 0;
    
//     while (true) {
        
//         for (int i = 0; i < n; i++) {
//             result[j][i] = vectors[i][indexes[i]];
//         }
//         j++;
        
//         size_t depth = 0;
//         while (depth < n && indexes[depth] >= sizes[depth]-1) {
//             depth++;
//         }
//         if (depth == n) {
//             break;
//         }
//         for (int i = 0; i < depth; i++) {
//             indexes[i] = 0;
//         }
//         indexes[depth] += 1;
        
//     }
    
//     return result;
// }


template <typename... Types, std::size_t... I>
vector<tuple<Types...>> cartesian_product_implementation(std::index_sequence<I...>,
                                                        const std::vector<Types>&... vec) {
    const static size_t n = sizeof...(Types);
    array<size_t, sizeof...(Types)> indexes {};
    array<size_t, sizeof...(Types)> sizes {vec.size()...};
    
    vector<tuple<Types...>> result;
    
    while (true) {
        result.push_back(make_tuple(vec[indexes[I]]...));
        
        size_t depth = 0;
        while (depth < n && indexes[depth] >= sizes[depth]-1) {
            depth++;
        }
        if (depth == n) {
            break;
        }
        for (size_t i = 0; i < depth; i++) {
            indexes[i] = 0;
        }
        indexes[depth] += 1;
    }   
    return result;
}

template <typename... Types>
vector<tuple<Types...>> cartesian_product(const std::vector<Types>&... vec) {
    return cartesian_product_implementation(std::index_sequence_for<Types...>{}, std::forward<decltype(vec)>(vec)...);
}




// Example usage
// int main() {
//     std::vector<std::vector<int>> vectors = {{1, 2}, {3, 4, 0}, {5}};
//     auto result = cartesian_product(vectors);

//     int ༼ = 0;
//     for (const auto& tuple : result) {
//         for (const auto& val : tuple) {
//             std::cout << val << " ";
//         }
//         std::cout << "\n";
//     }

//     return 0;
// }
