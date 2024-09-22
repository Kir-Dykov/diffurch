#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>
#include <thread>

#include "library/json.hpp"
#include "library/utils.hpp"
#include "library/save.hpp"
#include "library/discoque.hpp"
#include "library/dde.hpp"
#include "library/ddes_equations.hpp"
#include "library/runge_kutta_tables.hpp"
#include "library/progress_bar.hpp"

#include "library/vec.hpp"
#include "library/types.hpp"


using namespace std;

using json = nlohmann::json;

// helper template function
template <std::size_t args_n,  std::size_t... I>
auto unpack_json_doubles_sequence(json json_params, std::index_sequence<I...>, array<const char*, args_n>args) {
    return std::tuple((double)json_params[args[I]]...);
}
template <size_t args_n>
auto unpack_json_doubles(json json_params, array<const char*, args_n> args) {
    return unpack_json_doubles_sequence(json_params, std::make_index_sequence<args_n>{}, args);
}


int main(int argc, char* argv[]) {

    // feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);


    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	string params = argv[1];
	json json_params = json::parse(params);
	cout << "~~~  parameters: " << params << " ~~~" << endl;
    
    
	string output_filename = argv[2];
    
    
    auto [alpha_l, alpha_r, alpha_n, 
          tau_l, tau_r, tau_n, 
          t_finish, h, 
          A_l, A_r, A_n, 
          w_l, w_r, w_n] =
    unpack_json_doubles(json_params, array{"alpha_l", "alpha_r", "alpha_n",
                                           "tau_l", "tau_r", "tau_n", 
                                           "t_finish", "h",
                                           "A_l", "A_r", "A_n",
                                           "w_l", "w_r", "w_n"});
    
    vector<vector<double>> output(alpha_n, vector<double>(tau_n, 0));
    
    vector<double> alphas = linspace(alpha_l, alpha_r, alpha_n);
    vector<double> taus   = linspace(tau_l,     tau_r,   tau_n);
    vector<double> As   = linspace(A_l,     A_r,   A_n);
    vector<double> ws   = linspace(w_l,     w_r,   w_n);
    
    
    
    ProgressBar progress_bar(alpha_n*tau_n);
    
    
    vector<thread> threads;
    
    
    for (size_t alpha_i = 0; alpha_i < alpha_n; alpha_i++) {
        threads.push_back(thread([&, alpha_i](){
            for (size_t tau_i = 0; tau_i < tau_n; tau_i++) {
                double alpha = alphas[alpha_i];
                double tau   = taus[tau_i];

                auto de = ClipDDE1(alpha, tau);
                
                double discontinuity_count = 0.;
                
                auto event_function = [&discontinuity_count, t_finish](RK_TimeSeries<1, 1>& TS) {
                    thread_local double x_prev = 0; // thread_local is like static, but for each thread it has its own copy (with static the variable is shared across all threads)
                    double x = TS.x_ts.back();
                    if (TS.t_ts.back() > t_finish/2  // only count in the second half of a time series
                        && (x - 1.) * (x_prev - 1.) <= 0 // count zero crosses
                        && x != 1.) { // but not count them twice
                        discontinuity_count += 1.;
                    }
                    x_prev = x;
                };
                
                // VecFunc<1> phi            = [&](double t) {return Vec<1>{1.1};};
                // VecFunc<1> phi_derivative = [&](double t) {return Vec<1>{0};};
                
                VecFuncC<1,1> phi ([&](double t) {return Vec<1>{1.1};}, {[&](double t) {return Vec<1>{0};}});
                
                
                de.solution<ReturnOption::None,1>(h, t_finish, phi, event_function);
                
                vector<double> unique_counts = {discontinuity_count};
                
                double count_periodic_solutions = 0;
                if (discontinuity_count > 2) {
                    count_periodic_solutions+= 1.;
                }
                
                for (size_t A_i = 0; A_i < A_n; A_i++) {
                    for (size_t w_i = 0; w_i < w_n; w_i++) {
                        discontinuity_count = 0;
                        double A = As[A_i];
                        double w = ws[w_i];
                        
                        // VecFunc<1> phi            = 
                        //     [&](double t) {return Vec<1>{A*cos(w*t)};};
                        // VecFunc<1> phi_derivative = [&](double t) {return Vec<1>{-w*A*sin(w*t)};};
                        
                        VecFuncC<1,1> phi ([&](double t) {return Vec<1>{A*cos(w*t)};}, {[&](double t) {return Vec<1>{-w*A*sin(w*t)};}});
    
                        // output = de.solution<ReturnOption::AllInit>(h, t_finish, phi);
    
                        de.solution<ReturnOption::None,1>(h, t_finish, phi, event_function);
                        
                        bool is_unique = true;
                        for (auto C : unique_counts) {
                            if (abs(C - discontinuity_count) <= 2) {
                                is_unique = false;
                                break;
                            }
                        }
                        if (is_unique) {
                            // double checking for unstable solutions
                            
                            // count_periodic_solutions++;
                            // unique_counts.push_back(discontinuity_count);
                                
                            double prev_discontinuity_count = discontinuity_count;
                            discontinuity_count = 0;
                            de.solution<ReturnOption::None,1>(h, t_finish*2, phi, event_function); // we expect discontinuity_count be about twice the previous one if it is a stable cycle
                            
                            if (abs(prev_discontinuity_count - discontinuity_count/2) <= 4) {
                                count_periodic_solutions++;
                                unique_counts.push_back(prev_discontinuity_count);

                                //  if (count_periodic_solutions > 1) {
                                //     cout << count_periodic_solutions << ":    tau = " << tau << ",  alpha = " << alpha << ",    A = " << A << ",   w = " << w  << "," << endl;
                                // }
                            } else {
                                count_periodic_solutions+= 0.001;
                                // unique_counts.push_back(discontinuity_count);
                            }
                        
                            
                        }                      
                    }
                }
                
               
                

                output[alpha_i][tau_i] = count_periodic_solutions;         
                progress_bar.progress++;
                progress_bar.update();
                progress_bar.print();
            }
            
        })); 
        // threads.back().join();
    }
    
    for (int i = 0; i < threads.size(); i++) {
        threads[i].join();
    }
    
    
    
	save(output, "../output/bin/" + output_filename + ".bin");
    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;

    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
    return 0;
}