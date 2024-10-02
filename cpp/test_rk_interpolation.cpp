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

#include "library/real.hpp"

#include "equations/ode.hpp"
#include "equations/dde.hpp"

using namespace std;

// #ifndef DE
// #define DE ODE_lin_1
// #endif

template <int factor>
struct ReturnDenseSolution{
    template <size_t n, size_t phi_derivatives_n>
    static auto Return(RK_TimeSeries<n, phi_derivatives_n>& TS) { 
        vector<Real> t_dense_ts;
        vector<Vec<n>> x_dense_ts;
        Real t = 0;
        Real t_finish = TS.t_ts.back();
        Real dense_h = t_finish/TS.t_ts.size()/factor;
        
        while (t < t_finish) {
            t_dense_ts.push_back(t);
            x_dense_ts.push_back(TS.eval(t));
            t += dense_h;
        }
        
        return make_tuple(t_dense_ts, x_dense_ts);
    };
};

template <typename DE>
auto test_rk_interpolation(DE de, Real t_finish, vector<Real> hs) {
    int h_n = hs.size();
    
    vector<Real> error(h_n);
    
    for (int h_i = h_n-1; h_i >=0; h_i--) {
        Real h = hs[h_i];
        
        auto [t_dense_ts, x_dense_ts] = de.template solution<ReturnDenseSolution<100>>(h, t_finish, de.analytic_solutions[0]);
        auto true_x_dense_ts = de.analytic_solutions[0].eval_series(t_dense_ts);
        
        error[h_i] = 0;
        for (int i = 0; i < t_dense_ts.size(); i++) {
            error[h_i] = max(error[h_i], norm(x_dense_ts[i] - true_x_dense_ts[i]));
        } 
    }
    return error;
};



int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	string params = argv[1];
	json json_params = json::parse(params);
	cout << "~~~  parameters: " << params << " ~~~" << endl;
	string output_filename = argv[2];
    
    
    vector<vector<Real>> output;
    
    auto [t_finish, h] = json_unpack<Real, "t_finish", "h">(json_params);
    auto de = from_json<EQ>(json_params);   

    vector<Real> hs = expspace(0.01r, 1r, 100);
    auto er = test_rk_interpolation(de, t_finish, hs);
    
	save_arrays("../output/bin/" + output_filename + ".bin", hs, er);
    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;
    
    // auto x_last = de.template solution<???>(h, t_finish, de.analytic_solutions[0]);

    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}



