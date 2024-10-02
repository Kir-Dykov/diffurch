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

using namespace std;

// #ifndef DE
// #define DE ODE_lin_1
// #endif


template <typename DE>
tuple<vector<double>, vector<int>> test_rk(DE de, double t_finish, vector<double> hs) {
    int h_n = hs.size();
    vector<double> error(h_n);
    vector<int> time(h_n);
    for (int h_i = h_n-1; h_i >=0; h_i--) {
        double h = hs[h_i];
        
        auto clock_begin = std::chrono::high_resolution_clock::now();
        auto x_last = de.template solution<ReturnLastValue>(h, t_finish, de.analytic_solutions[0]);
        // auto x_last = de.template solution<ReturnOption::Last>(h, t_finish, de.analytic_solutions[0]);
        auto x_last_true = de.analytic_solutions[0](t_finish);
        auto diff = (x_last - x_last_true);
        auto clock_end = std::chrono::high_resolution_clock::now();
        
        int nanoseconds = chrono::duration_cast<chrono::nanoseconds>(clock_end - clock_begin).count();
        
        error[h_i] = sqrt(diff*diff);
        time[h_i] = nanoseconds;
    }
    return make_tuple(error, time);
}



int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	string params = argv[1];
	json json_params = json::parse(params);
	cout << "~~~  parameters: " << params << " ~~~" << endl;
	string output_filename = argv[2];
    
    
    vector<vector<double>> output;
    
    auto [t_finish, h] = json_unpack<double, "t_finish", "h">(json_params);
    auto de = from_json<EQ>(json_params);   

    vector<double> hs = expspace(0.01, 1., 1000);
    auto [er, time] = test_rk(de, t_finish, hs);
    
	save_arrays("../output/bin/" + output_filename + ".bin", hs, er, time);
    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;
    
    // auto x_last = de.template solution<???>(h, t_finish, de.analytic_solutions[0]);

    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}



