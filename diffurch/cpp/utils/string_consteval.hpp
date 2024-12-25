// #pragma once

// #include <algorithm>

// // code were taken from https://dev.to/sgf4/strings-as-template-parameters-c20-4joh

// template<std::size_t N>
// struct string_consteval {
//     char data[N] {};

//     consteval string_consteval(const char (&str)[N]) {
//         std::copy_n(str, N, data);
//     }

//     consteval bool operator==(const string_consteval<N> str) const {
//         return std::equal(str.data, str.data+N, data);
//     }

//     template<std::size_t N2>
//     consteval bool operator==(const string_consteval<N2> s) const {
//         return false;
//     }

//     template<std::size_t N2>
//     consteval string_consteval<N+N2-1> operator+(const string_consteval<N2> str) const {
//         char newchar[N+N2-1] {};
//         std::copy_n(data, N-1, newchar);
//         std::copy_n(str.data, N2, newchar+N-1);
//         return newchar;
//     }

//     consteval char operator[](std::size_t n) const {
//         return data[n];
//     }

//     consteval std::size_t size() const {
//         return N-1;
//     }
// };

// template<std::size_t s1, std::size_t s2>
// consteval auto operator+(string_consteval<s1> fs, const char (&str) [s2]) {
//     return fs + string_consteval<s2>(str);
// }

// template<std::size_t s1, std::size_t s2>
// consteval auto operator+(const char (&str) [s2], string_consteval<s1> fs) {
//     return string_consteval<s2>(str) + fs;
// }

// template<std::size_t s1, std::size_t s2>
// consteval auto operator==(string_consteval<s1> fs, const char (&str) [s2]) {
//     return fs == string_consteval<s2>(str);
// }

// template<std::size_t s1, std::size_t s2>
// consteval auto operator==(const char (&str) [s2], string_consteval<s1> fs) {
//     return string_consteval<s2>(str) == fs;
// }