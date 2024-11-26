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


#include "equations/ode.hpp"
#include "equations/dde.hpp"
#include "equations/relay_dde.hpp"

using namespace std;

// #ifndef DE
// #define DE ODE_lin_1
// #endif



int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	// cout << "~~~  parameters: " << params << " ~~~" << endl;    

    json JSON_PARAMS = json::parse(argv[1]);
	string output_filename = argv[2];
    
    JSON_UNPACK(Real, t_finish, h);    
    auto de = from_json<EQ>(JSON_PARAMS);   
    
    auto [t, x] = de.solution(h, t_finish, de.analytic_solutions[0], ReturnDenseSolution(100));
    auto true_x = de.analytic_solutions[0].eval_series(t);
    save_arrays("../output/bin/" + output_filename + ".bin", t, x, true_x);

    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;
    
    // auto x_last = de.template solution<???>(h, t_finish, de.analytic_solutions[0]);

    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}



