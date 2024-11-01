#include <iostream>

#include "../library/include.hpp"
#include "dde_clip_2_equation.hpp"
#include "../library/switches.hpp"

using namespace std;

/*
for each initial condition of the form A*cos(w*t)
we want to 
* compute the period of the solution, to do that
    ** make sequence of switch-packs
    ** starting from the second to last, we seach for the matching switch-pack
    ** period is the difference between time of matched switch-packs
* then, if such period has not been found, we add it to the
*/


int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	string params = argv[1];
	json json_params = json::parse(params);
	cout << "~~~  parameters: " << params << " ~~~" << endl;
	string output_filename = argv[2];
    
    
    auto [t_finish, h] = json_unpack<Real, "t_finish", "h">(json_params);
    auto [A_l, A_r, w_l, w_r] = json_unpack<Real, "A_l", "A_r", "w_l", "w_r">(json_params);
    auto [A_n, w_n] = json_unpack<int, "A_n", "w_n">(json_params);
    // auto [T] = json_unpack<Real, "T">(json_params);
    auto [tau] = json_unpack<Real, "tau">(json_params);
    
    
    auto de = from_json<DDE_clip_2>(json_params);   
    
    vector<vector<Real>> output(A_n, vector<Real>(w_n, 0));
    
    // vector<double> alphas = linspace(alpha_l, alpha_r, alpha_n);
    // vector<double> taus   = linspace(tau_l,     tau_r,   tau_n);
    vector<double> As   = linspace(A_l,     A_r,   A_n);
    vector<double> ws   = linspace(w_l,     w_r,   w_n);
    
    auto return_period = ReturnPeriod<2>([](Vec<2> X) -> Real {return X[0];}, {-1.r, 1.r}, tau);
    
    for (int i = 0; i < A_n; i++) {
        for (int j = 0; j < w_n; j++) {
            auto A = As[i];
            auto w = ws[j];
            
            Real period = de.solution(h, t_finish, Cos2VecMapC(A, w), return_period);
            output[i][j] = period;
        }
    }
        
	save_arrays("../output/bin/" + output_filename + ".bin", As, ws, output);
    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;
    
    // auto x_last = de.template solution<???>(h, t_finish, de.analytic_solutions[0]);

    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}



