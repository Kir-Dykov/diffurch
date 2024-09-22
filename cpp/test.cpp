#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>
#include <chrono>
#include <tuple>
#include <thread>

#include "library/json.hpp"

#include "library/vec.hpp"
// #include "library/dde.hpp"
// #include "library/ddes_equations.hpp"

#include "library/types.hpp"

using namespace std;

template<typename T>
constexpr bool is_default(T x) {
    return x == T{};
}


size_t verbose_square(size_t x) {
    cout << "Squaring " << x << endl;
    return x*x;
}


template <size_t degree>
struct test {
    
    template <size_t order = 0>
    void eval() {
        static size_t square = verbose_square(order);
        
        square += 1000;
        cout << square << endl;
    }
};

template <size_t n>
struct array_test_class {
    array<double, n> coefs;
    
    array_test_class(array<double, n> coefs) : coefs(coefs) {};
};




// void thread_lambda() {
//     auto count_lambda = []() {
//         static int count = 0;  // Static variable
//         count++;               // Increment
//         std::cout << "Thread ID: " << std::this_thread::get_id() << " | Count: " << count << std::endl;
//     };

//     count_lambda();
// }


int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    
    
    array_test_class qwe {{1., 2., 3.}};
    
    test<1> x;
    
    x.eval();
    x.eval<1>();
    x.eval<2>();
    x.eval<1>();
    
    test<1> y;
    
    y.eval();
    y.eval<1>();
    
    // VecFunc<3> x{};
    // cout << *x << endl;
    // cout << (void*)x.target() << endl;
    // if constexpr (x) {
    //     cout << "HUA" << endl;
    // } else {
    //     cout << "HUI" << endl;
    // }
    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}



