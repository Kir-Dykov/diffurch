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
    // auto [A, w, gamma] = json_unpack<Real, "A", "w", "gamma">(json_params);
    
    auto [T, alpha] = json_unpack<Real, "T", "alpha">(json_params);
    
    
    auto de = from_json<EQ>(json_params);   
    auto de_relay = from_json<RelayNDDE1>(json_params);   
    
    
   
    
    auto [k] = json_unpack<int, "k">(json_params);
        
    // vector<Real> t_switches = {0., T};
    vector<vector<Real>> t_switches_options;
    t_switches_options.push_back({0., T}); // slow oscillations 0 extra switches
    t_switches_options.push_back({0., T/8, T}); // + 1
    t_switches_options.push_back({0., T/8, 7*T/8, T}); // + 1
    t_switches_options.push_back({0.,T/6, T/5, T/4, T/3, T/2, 2*T/3, 3*T/4, 4*T/5, 5*T/6, T});
    
    vector<Real> t_switches = t_switches_options[k];
    auto phi = de_relay.analytic_periodic_solution(t_switches);
    
    
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
    

    
    

    // de.eps = eps0;
     #ifdef USE_SMOOTH_SOLUTION
    auto [cos_a, cos_w, cos_gamma] = json_unpack<Real, "cos_a", "cos_w", "cos_gamma">(json_params);
    auto phi_cos = CosVecMapC<1>(cos_a, cos_w, cos_gamma);
    auto [ts, xs, dot_xs] = de.template solution<ReturnSolutionAndDerivative>(h, t_finish, phi_cos);
    #else
    auto [ts, xs, dot_xs] = de.template solution<ReturnSolutionAndDerivative>(h, t_finish, phi);
    #endif
    
    auto reference = phi.eval_series(ts);
    auto dot_reference = phi.eval_series<1>(ts);
    
    vector<Real> zs = linspace(-5.r, 5.r, 10000);
    vector<Real> fs (zs.size());
    std::transform(zs.begin(), zs.end(), fs.begin(), [&de](Real z){return de.F_eps(z);});

    #ifdef USE_SMOOTH_SOLUTION
    vector<Real> init_t = linspace(-T, 0r, 1000);
    vector<Real> init_x (init_t.size());
    vector<Real> init_dot_x (init_t.size());
    std::transform(init_t.begin(), init_t.end(), init_x.begin(),     [&phi_cos](Real t){return phi_cos.eval(t);});
    std::transform(init_t.begin(), init_t.end(), init_dot_x.begin(), [&phi_cos](Real t){return phi_cos.template eval<1>(t);});
    
	save_arrays("../output/bin/" + output_filename + ".bin", ts, xs, dot_xs, reference, dot_reference, zs, fs, init_t, init_x, init_dot_x);
    #else
	save_arrays("../output/bin/" + output_filename + ".bin", ts, xs, dot_xs, reference, dot_reference, zs, fs);
    #endif
    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;
    
    // auto x_last = de.template solution<???>(h, t_finish, de.analytic_solutions[0]);

    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}



