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

#include "library/json_unpack.hpp"
#include "library/json.hpp"
// using json = nlohmann::json;
#include "library/save.hpp"
#include "library/vec.hpp"


// #include <boost/hana.hpp>
// #include <type_traits>
// namespace hana = boost::hana;

// #include "equations/ode.hpp"
// #include "equations/dde.hpp"

using namespace std;


// #define DEFAULT "default"


// #include "library/string_consteval.hpp"
// #include "library/param_names.hpp"
// #include "library/tuple.hpp"

// #include "library/macros.hpp"
// template <string_consteval s, typename P>
// using param_index_v = param_index<s, P>::value;

#include "library/real.hpp"

// #ifndef X
// #define X alpha
// #endif

// template <typename T, typename U>
// struct tuple_concat {};

// template <typename... Args1, typename... Args2>
// struct tuple_concat<tuple<Args1...>, tuple<Args2...>> {
//     using type = tuple<Args1..., Args2...>;
// };



// template <typename T, typename U>
// using tuple_concat_t = tuple_concat<T, U>::type;



// template <size_t s>
// struct RK_table {
//     using A_type = tuple_concat_t< typename RK_table<s-1>::A_type, tuple<array<Real, s-1>> >;
//     using B_type = array<Real, s>;
//     using C_type = array<Real, s>;
    
// };

// template <>
// struct RK_table<0> {
//     using A_type = tuple<>;
// };



// struct Euler : RK_table<1> {
//     A_type A = {{}};
//     B_type B = {1.};
//     C_type C = {0.};
// };

// struct RK4 : RK_table<4> {
//     A_type A = 
//         {{}, 
//          {1./2.}, 
//          {0.,     1./2.}, 
//          {0.,     0.,      1.}
//         };
//     B_type B = {1./6., 1./3., 1./3., 1./6.};
//     C_type C = 
//     {
//         0,
//         1./2.,
//         1./2.,
//         1.
//     };
// };

#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <quadmath.h>

using namespace boost::multiprecision;




int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	json json_params = json::parse(argv[1]);
	cout << "~~~  parameters: " << argv[1] << " ~~~" << endl;
	string output_filename = argv[2];
    
    
    // using boost::multiprecision::float128;
    
    
    auto n = 123r;
    cout << typeid(decltype(n)).name() << endl;
    
    cout << sizeof(double) << endl;
    cout << sizeof(Real) << endl;
//     cout << q1 << endl;
//     cout << q2 << endl;
//     cout << q3 << endl;

    // using t = tuple_concat_t<tuple<>, tuple<int>>;
    
    
    // RK4 rk;
    // // cout << get<3>(rk.A) << endl;
    // for (int i = 0; i < 4; i++) {
    //     cout << get<i>(rk.A) << endl;
    // }

    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;
    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}



