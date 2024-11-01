#include <iostream>

#include "../library/include.hpp"
#include "dde_clip_2_equation.hpp"


using namespace std;


// template <typename DE>
// auto test_rk_interpolation(DE de, Real t_finish, vector<Real> hs) {
//     int h_n = hs.size();
    
//     vector<Real> error(h_n);
    
//     for (int h_i = h_n-1; h_i >=0; h_i--) {
//         Real h = hs[h_i];
        
//         auto [t_dense_ts, x_dense_ts] = de.template solution<ReturnDenseSolution<100>>(h, t_finish, de.analytic_solutions[0]);
//         auto true_x_dense_ts = de.analytic_solutions[0].eval_series(t_dense_ts);
        
//         error[h_i] = 0;
//         for (int i = 0; i < t_dense_ts.size(); i++) {
//             error[h_i] = max(error[h_i], norm(x_dense_ts[i] - true_x_dense_ts[i]));
//         } 
//     }
//     return error;
// };



int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	string params = argv[1];
	json json_params = json::parse(params);
	cout << "~~~  parameters: " << params << " ~~~" << endl;
	string output_filename = argv[2];
    
    
    vector<vector<Real>> output;
    
    auto [t_finish, h] = json_unpack<Real, "t_finish", "h">(json_params);
    auto [A, w] = json_unpack<Real, "A", "w">(json_params);
    auto [tau] = json_unpack<Real, "tau">(json_params);
    
    auto de = from_json<DDE_clip_2>(json_params);   
    
    // cout << (has_Return<ReturnSolution, RK_TimeSeries<1, 0>>) << endl;
    
    auto [ts, xs] = de.solution(h, t_finish, Cos2VecMapC(A, w), ReturnSolution{});
    
    Real period = de.solution(h, t_finish, Cos2VecMapC(A, w), ReturnPeriod<2>([](Vec<2> X) -> Real {return X[0];}, {-1.r, 1.r}, tau));
    
    auto [ts_p, xs_p] = de.solution(h, t_finish, Cos2VecMapC(A, w), ReturnPoincareSequence<2>([](Vec<2> X) -> Real {return X[0];}, {-1.r, 1.r}, tau));
    
    cout << "period = " << period << endl;
    
    vector<Real> pp = {period};
    
    save_arrays("../output/bin/" + output_filename + ".bin", ts, xs, ts_p, xs_p, pp);
    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;
    
    // auto x_last = de.template solution<???>(h, t_finish, de.analytic_solutions[0]);

    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}



