#include <iostream>

#include "../library/include.hpp"
#include "ndde_relay_equations.hpp"


using namespace std;




int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	string params = argv[1];
	json json_params = json::parse(params);
	cout << "~~~  parameters: " << params << " ~~~" << endl;
	string output_filename = argv[2];
    
    
    auto [t_finish, h] = json_unpack<Real, "t_finish", "h">(json_params);
    auto [A, w, gamma] = json_unpack<Real, "A", "w", "gamma">(json_params);
    
    auto [T, alpha] = json_unpack<Real, "T", "alpha">(json_params);
    
    
    // auto de = from_json<RelayNDDE1arctan>(json_params);   
    auto de = from_json<RelayNDDE1>(json_params);   
    
    // auto [ts, xs] = de.template solution<ReturnSolution>(h, t_finish, CosVecMapC<1>(A, w, gamma));
    vector<Real> t_switches = {0.,T/6, T/5, T/4, T/3, T/2, 2*T/3, 3*T/4, 4*T/5, 5*T/6, T};
    
    auto phi = de.analytic_periodic_solution(t_switches);
    auto [ts, xs] = de.template solution<ReturnSolution>(h, t_finish, phi);
    auto [ts2, xs2] = de.template solution<ReturnSolution>(h, t_finish, VecMapC<1,1,1>([&phi, A](Real t){return A + phi.eval(t);}, phi.df));
    
    auto true_xs = phi.eval_series(ts);
    
    auto     X_points = phi.eval_series<0>(t_switches);
    auto dot_X_points = phi.eval_series<1>(t_switches);
    Real max_x = std::transform_reduce(X_points.begin(), X_points.end(), 0.r, [](Real a, Real b) { return std::max(a, b); }, [](Real a){ return std::abs(a);});
    Real min_dot_x = std::transform_reduce(dot_X_points.begin(), dot_X_points.end(), 999999999.r, [](Real a, Real b) { return std::min(a, b); }, [](Real a){ return std::abs(a);});
    
    cout << X_points << endl;
    cout << max_x << endl;
     cout << dot_X_points << endl;
    cout << min_dot_x << endl;
    
    Real delta = min(min_dot_x/abs(alpha), 0.5*(1 - max_x/abs(alpha)));
    Real eps0 = abs(alpha)*delta/tan((1-delta*0.5)*M_PI*0.5);
    cout << "delta = " << delta << endl;
    cout << "epsilon = " << eps0 << endl;



	save_arrays("../output/bin/" + output_filename + ".bin", ts, xs, ts2, xs2, true_xs);
    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;
    
    // auto x_last = de.template solution<???>(h, t_finish, de.analytic_solutions[0]);

    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}



