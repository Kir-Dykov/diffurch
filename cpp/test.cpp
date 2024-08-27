#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>

#include "json.hpp"

#include "array_operations.hpp"

using namespace std;


// template<int N, typename Enable = void>
// struct Vec;
// template<int N>
// struct Vec<N, typename std::enable_if<N == 1>::type> {
//     using Type = double;
// };
// template<int N>
// struct Vec<N, typename std::enable_if<(N > 1)>::type> {
//     using Type = array<double, N>;
// };
// template<int N>
// struct Vec<N, typename std::enable_if<(N < 1)>::type> {
//     using Type = void;
// };

// template <int N, int M>
// using Func = function<Vec<M>(const Vec<N>&)>;
// template <int N1, int N2, int M>
// using Func2 = function<Vec<M>(const Vec<N1>&, const Vec<N2>&)>;


int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;

    Vec<1> x = 2.;
    
    cout << x << endl;
    
    Vec<2> y = {3., 4.};
    
    cout << y << endl;
    
    Func<1,1> half = [](double x){return .5*x;};
    
    Func<2,1> sum = [](array<double, 2> x){return x[0] + x[1];};
    
    cout << half(x) << endl;
    
    cout << sum(y) << endl;
    
    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}
