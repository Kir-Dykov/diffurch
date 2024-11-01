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


#include "../library/json_unpack.hpp"
#include "../library/json.hpp"

using json = nlohmann::json;
// #include "../library/save.hpp"
// #include "../library/vec.hpp"

// #include "equations/ode.hpp"
// #include "equations/dde.hpp"

using namespace std;


// #include "../library/real.hpp"

// #include <math.h>




template <typename T>
concept has_param_names = requires {
    typename T::param_names;
};

class MyClass1 {
    using param_names = int;
}; 

class MyClass0 {
    
};


template <typename T>
concept HasParamNames = requires {
    typename T::param_names; // Checks if T has a nested type 'param_names'
};

// Example class that satisfies the concept
class MyClass {
public:
    using param_names = int; // Defines the type alias 'param_names'
};

// Example class that does NOT satisfy the concept
class AnotherClass {
    // No 'param_names' type alias defined
};

// Test function that accepts only classes satisfying HasParamNames concept
template <HasParamNames T>
void checkParamNames(T) {
    // Do something knowing that T has 'param_names'
}





int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	json json_params = json::parse(argv[1]);
	cout << "~~~  parameters: " << argv[1] << " ~~~" << endl;
	string output_filename = argv[2];
    
        MyClass obj1;
    checkParamNames(obj1); // This will compile
    
    cout << has_param_names<MyClass0> << endl;
    cout << has_param_names<MyClass1> << endl;

    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;
    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}



