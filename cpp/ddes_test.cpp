#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>

#include <cfenv>

// #include <boost/numeric/interval.hpp>
// using namespace boost::numeric;

#include "library/json.hpp"
#include "library/utils.hpp"
#include "library/save.hpp"
#include "library/discoque.hpp"
#include "library/dde.hpp"
#include "library/ddes_equations.hpp"
#include "library/runge_kutta_tables.hpp"

#include "library/vec.hpp"
#include "library/types.hpp"


using namespace std;

using json = nlohmann::json;

double b,c,d,tau,v,T,h;
string eq;


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
	string output_filename = argv[2];
    
    // cout << "~~~ output_filename : " << output_filename << " ~~~" << endl;
    
	json json_params = json::parse(params);

	cout << "~~~  parameters: " << params << " ~~~" << endl;
    
    if (!json_params.contains("eq")) {
        cout << "No 'eq' parameter specifyed, exiting." << endl;
        return 1;
    }
    
    
    
    eq = (string)json_params["eq"];
    
    
    vector<vector<double>> output;
    
        

    if (eq == "Relay2") {        
        auto [b, c, d, tau, v0, x0, t_finish, h] =
        unpack_json_doubles(json_params, array{"b", "c", "d", "tau", "v0", "x0", "t_finish", "h"});
        
        auto de = SignDDE2(b, c, d, tau);
        
        VecFunc<2> phi            = [&](double t) {return Vec<2>{x0 + v0*t, v0};};
        VecFunc<2> phi_derivative = [&](double t) {return Vec<2>{v0, 0};};
        
        // order test:
//         if (t_finish < 20) {
//             double H = 0.5;
//             for (int i = 0; i < 10; i++) {
//                 auto x0 = de.solution<ReturnOption::Last>(H,   t_finish, phi, phi_derivative)[0];
//                 auto x1 = de.solution<ReturnOption::Last>(H/2, t_finish, phi, phi_derivative)[0];
//                 auto x2 = de.solution<ReturnOption::Last>(H/4, t_finish, phi, phi_derivative)[0];
                
//                 cout << "H = " << H << " and error ratio = " << (x1-x0)/(H*H*H*H) << endl;
                
//                 H = H/2;
//             }
//         }
        
        output = de.solution<ReturnOption::All>(h, t_finish, phi, phi_derivative);
    } 
    else if (eq == "Relay1") {
        auto [c, d0, d1, tau, t_finish, h] =
        unpack_json_doubles(json_params, array{"c", "d0", "d1", "tau", "t_finish", "h"});
       
        auto de = SignDDE1(c, d0, d1, tau);
        
        VecFunc<1> phi            = [](double t) {return Vec<1>{t};};
        VecFunc<1> phi_derivative = [](double t) {return Vec<1>{1};};
        
        
        
        // order test:
        
//         if (t_finish < 20) {
//             double H = 0.5;
//             for (int i = 0; i < 10; i++) {
//                 auto x0 = de.solution<ReturnOption::Last>(H,   t_finish, phi, phi_derivative)[0];
//                 auto x1 = de.solution<ReturnOption::Last>(H/2, t_finish, phi, phi_derivative)[0];
//                 auto x2 = de.solution<ReturnOption::Last>(H/4, t_finish, phi, phi_derivative)[0];
                
//                 cout << "H = " << H << " and error ratio = " << (x1-x0)/(H*H) << endl;
                
//                 H = H/2;
//             }
//         }
        
        
        output = de.solution<ReturnOption::All>(h, t_finish, phi, phi_derivative);
        
    }  
//     else if (eq == "Croissant") {
//         // x' + b x = Sign( x_tau )
//         double b = (double)json_params["b"];
//         double d0= (double)json_params["d0"];
//         double d1= (double)json_params["d1"];
//         double tau=(double)json_params["tau"];
//         double t_finish = (double)json_params["t_finish"];
//         double h = (double)json_params["h"];
        
//         auto dde = Croissant(b, d0, d1, tau);
//         Func<1,1> phi            = [](double t) {return t;};
//         Func<1,1> phi_derivative = [](double t) {return 1;};
//         InitialValueProblem<1,1,1> ivp (dde, phi, phi_derivative);
//         // ivp.compute_solution(h, t_finish);
//         Func<1,1> psi = [](double t){return 1.;};
//         ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
//         output = ivp.get_solution_and_lyapunov_exponent();     
//     }
//     else if (eq == "MakeyGlass") {
//         // x' + b x = Sign( x_tau )
//         double b = (double)json_params["b"];
//         double d0= (double)json_params["d0"];
//         double d1= (double)json_params["d1"];
//         double tau=(double)json_params["tau"];
//         double x0 = (double)json_params["x0"];
//         double v0 = (double)json_params["v0"];
        
        
//         double t_finish = (double)json_params["t_finish"];
//         double h = (double)json_params["h"];
        
//         auto dde = MakeyGlass(b, d0, d1, tau);
//         Func<1,1> phi            = [x0,v0,tau](double t) {return exp(x0*cos(2*M_PI*v0/ tau * t ));};
//         Func<1,1> phi_derivative = [x0,v0,tau](double t) {return -x0*2*M_PI*v0/tau*sin(2*M_PI*v0/ tau * t )*exp(x0*cos(2*M_PI*v0/ tau * t ));};
//         InitialValueProblem<1,1,1> ivp (dde, phi, phi_derivative);
//         // ivp.compute_solution(h, t_finish);
//         Func<1,1> psi = [](double t){return 1.;};
//         ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
//         output = ivp.get_solution_and_lyapunov_exponent();     
//     }
//     else if (eq == "MakeyGlassExp") {
//         // x' + b x = Sign( x_tau )
//         double b = (double)json_params["b"];
//         double d0= (double)json_params["d0"];
//         double d1= (double)json_params["d1"];
//         double tau=(double)json_params["tau"];
//         double x0 = (double)json_params["x0"];
//         double v0 = (double)json_params["v0"];
        
        
//         double t_finish = (double)json_params["t_finish"];
//         double h = (double)json_params["h"];
        
//         auto dde = MakeyGlassExp(b, d0, d1, tau);
//         Func<1,1> phi            = [x0,v0,tau](double t) {return x0*cos(2*M_PI*v0/ tau * t );};
//         Func<1,1> phi_derivative = [x0,v0,tau](double t) {return -x0*2*M_PI*v0/tau*sin(2*M_PI*v0/ tau * t );};
//         InitialValueProblem<1,1,1> ivp (dde, phi, phi_derivative);
//         // ivp.compute_solution(h, t_finish);
//         Func<1,1> psi = [](double t){return 1.;};
//         ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
//         output = ivp.get_solution_and_lyapunov_exponent();     
//     }
//      else if (eq == "Relay2damp") {
//         // x' + b x = Sign( x_tau )
//         double c = (double)json_params["c"];
//         double gamma= (double)json_params["gamma"];
//         double b0= (double)json_params["b0"];
//         double b1= (double)json_params["b1"];
//         double tau=(double)json_params["tau"];
//         double t_finish = (double)json_params["t_finish"];
//         double h = (double)json_params["h"];
        
//         auto dde = Relay2damp(gamma, b0, b1, c, tau);
//         Func<1,2> phi            = [&](double t) {return Vec<2>{0.5*d/c, 0};};
//         Func<1,2> phi_derivative = [&](double t) {return Vec<2>{0, 0};};
//         InitialValueProblem<2,1,0> ivp (dde, phi, phi_derivative);
//         // ivp.compute_solution(h, t_finish);
//         Func<1,2> psi = [tau](double t){return Vec<2>{sin(1.1*2*M_PI/tau*t), 1.1*2*M_PI/tau*cos(2*M_PI/tau*t)};};
//         ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
//         output = ivp.get_solution_and_lyapunov_exponent();     
//     }
//     else  if (eq == "RandomWalk") {
//         // x' + b x = Sign( x_tau )
//         double b = (double)json_params["b"];
//         double d0= (double)json_params["d0"];
//         double d1= (double)json_params["d1"];
//         double tau=(double)json_params["tau"];
//         double x0 = (double)json_params["x0"];
//         double t_finish = (double)json_params["t_finish"];
//         double h = (double)json_params["h"];
        
//         auto dde = RandomWalk(b, d0, d1, tau);
//         Func<1,1> phi            = [x0](double t) {return x0;};
//         Func<1,1> phi_derivative = [](double t) {return 0;};
//         InitialValueProblem<1,1,0> ivp (dde, phi, phi_derivative);
//         // ivp.compute_solution(h, t_finish);
//         Func<1,1> psi = [](double t){return 1.;};
//         ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
//         output = ivp.get_solution_and_lyapunov_exponent();
//     } 
    else {
        throw("Selected equation is not defined.");
    }
    
    
	save(output, "../output/bin/" + output_filename + ".bin");
    
	// save(output, "solution.bin");
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;

    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}

        // output = vector<vector<double>>{t_ts, x_ts};
        
        // InitialValueProblem<2,1,0> ivp (dde, phi, phi_derivative);
        
        // VecFunc<2> psi = [tau](double t){return Vec<2>{sin(1.1*2*M_PI/tau*t), 1.1*2*M_PI/tau*cos(2*M_PI/tau*t)};};
        // ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
        // output = ivp.get_solution_and_lyapunov_exponent();      