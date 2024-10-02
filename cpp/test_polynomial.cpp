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
using json = nlohmann::json;
#include "library/save.hpp"
#include "library/param_names.hpp"


#include "library/utils.hpp"

#include "library/polynomial.hpp"
// #include "equations/ode.hpp"
// #include "equations/dde.hpp"

using namespace std;

// #ifndef DE
// #define DE ODE_lin_1
// #endif




int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	string params = argv[1];
	json json_params = json::parse(params);
	cout << "~~~  parameters: " << params << " ~~~" << endl;
	string output_filename = argv[2];
    
    
    auto [a0, a1, a2, a3] = json_unpack<Real, "a0", "a1", "a2", "a3">(json_params);
    Polynomial<3> p(array<Real,4>{a0,a1,a2,a3});
    
    auto ts = linspace(0r, 1r, 1000);
    vector<Real> xs;
    for (int i = 0; i < ts.size(); i++) {
        xs.push_back(p.eval<0>(ts[i]));
    }

  
    
	save_arrays("../output/bin/" + output_filename + ".bin", ts, xs);
    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;
    
    // auto x_last = de.template solution<???>(h, t_finish, de.analytic_solutions[0]);

    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}



