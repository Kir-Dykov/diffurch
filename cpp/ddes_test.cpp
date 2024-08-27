#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>

#include "json.hpp"
#include <boost/numeric/interval.hpp>
using namespace boost::numeric;

#include "utils.hpp"
#include "save.hpp"
#include "discoque.hpp"
// #include "dde1.hpp"
#include "ddes.hpp"
#include "ddes_equations.hpp"
#include "array_operations.hpp"


using namespace std;

using json = nlohmann::json;

double b,c,d,tau,v,T,h;
string eq;


int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


	string params = argv[1];
	string output_prefix = argv[2];
	json json_params = json::parse(params);

	cout << "~~~  parameters: " << params << " ~~~" << endl;
    
    if (!json_params.contains("eq")) {
        cout << "No 'eq' parameter specifyed, exiting." << endl;
        return 1;
    }
    
    
    

	
    
    eq = (string)json_params["eq"];
    
    
    vector<vector<double>> output;
    
        
    if (eq == "Relay1") {
        // x' + b x = Sign( x_tau )
        double b = (double)json_params["b"];
        double d0= (double)json_params["d0"];
        double d1= (double)json_params["d1"];
        double tau=(double)json_params["tau"];
        double t_finish = (double)json_params["t_finish"];
        double h = (double)json_params["h"];
        
        auto dde = Relay1(b, d0, d1, tau);
        Func<1,1> phi            = [](double t) {return t;};
        Func<1,1> phi_derivative = [](double t) {return 1;};
        InitialValueProblem<1,1,0> ivp (dde, phi, phi_derivative);
        // ivp.compute_solution(h, t_finish);
        Func<1,1> psi = [](double t){return 1.;};
        ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
        output = ivp.get_solution_and_lyapunov_exponent();
        
    } else if (eq == "Relay2") {
        double b = (double)json_params["b"];
        double c = (double)json_params["c"];
        double d = (double)json_params["d"];
        double tau=(double)json_params["tau"];
        double v0 = (double)json_params["v0"];
        double x0 = (double)json_params["x0"];
        double t_finish = (double)json_params["t_finish"];
        double h = (double)json_params["h"];
        
        auto dde = Relay2(b, c, d, tau);
        
        Func<1,2> phi            = [&](double t) {return Vec<2>{x0 + v0*t, v0};};
        Func<1,2> phi_derivative = [&](double t) {return Vec<2>{v0, 0};};
        
        InitialValueProblem<2,1,0> ivp (dde, phi, phi_derivative);
        
        Func<1,2> psi = [tau](double t){return Vec<2>{sin(1.1*2*M_PI/tau*t), 1.1*2*M_PI/tau*cos(2*M_PI/tau*t)};};
        ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
        output = ivp.get_solution_and_lyapunov_exponent();        
    }
    else if (eq == "Croissant") {
        // x' + b x = Sign( x_tau )
        double b = (double)json_params["b"];
        double d0= (double)json_params["d0"];
        double d1= (double)json_params["d1"];
        double tau=(double)json_params["tau"];
        double t_finish = (double)json_params["t_finish"];
        double h = (double)json_params["h"];
        
        auto dde = Croissant(b, d0, d1, tau);
        Func<1,1> phi            = [](double t) {return t;};
        Func<1,1> phi_derivative = [](double t) {return 1;};
        InitialValueProblem<1,1,1> ivp (dde, phi, phi_derivative);
        // ivp.compute_solution(h, t_finish);
        Func<1,1> psi = [](double t){return 1.;};
        ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
        output = ivp.get_solution_and_lyapunov_exponent();     
    }
    else if (eq == "MakeyGlass") {
        // x' + b x = Sign( x_tau )
        double b = (double)json_params["b"];
        double d0= (double)json_params["d0"];
        double d1= (double)json_params["d1"];
        double tau=(double)json_params["tau"];
        double x0 = (double)json_params["x0"];
        double v0 = (double)json_params["v0"];
        
        
        double t_finish = (double)json_params["t_finish"];
        double h = (double)json_params["h"];
        
        auto dde = MakeyGlass(b, d0, d1, tau);
        Func<1,1> phi            = [x0,v0,tau](double t) {return exp(x0*cos(2*M_PI*v0/ tau * t ));};
        Func<1,1> phi_derivative = [x0,v0,tau](double t) {return -x0*2*M_PI*v0/tau*sin(2*M_PI*v0/ tau * t )*exp(x0*cos(2*M_PI*v0/ tau * t ));};
        InitialValueProblem<1,1,1> ivp (dde, phi, phi_derivative);
        // ivp.compute_solution(h, t_finish);
        Func<1,1> psi = [](double t){return 1.;};
        ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
        output = ivp.get_solution_and_lyapunov_exponent();     
    }
    else if (eq == "MakeyGlassExp") {
        // x' + b x = Sign( x_tau )
        double b = (double)json_params["b"];
        double d0= (double)json_params["d0"];
        double d1= (double)json_params["d1"];
        double tau=(double)json_params["tau"];
        double x0 = (double)json_params["x0"];
        double v0 = (double)json_params["v0"];
        
        
        double t_finish = (double)json_params["t_finish"];
        double h = (double)json_params["h"];
        
        auto dde = MakeyGlassExp(b, d0, d1, tau);
        Func<1,1> phi            = [x0,v0,tau](double t) {return x0*cos(2*M_PI*v0/ tau * t );};
        Func<1,1> phi_derivative = [x0,v0,tau](double t) {return -x0*2*M_PI*v0/tau*sin(2*M_PI*v0/ tau * t );};
        InitialValueProblem<1,1,1> ivp (dde, phi, phi_derivative);
        // ivp.compute_solution(h, t_finish);
        Func<1,1> psi = [](double t){return 1.;};
        ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
        output = ivp.get_solution_and_lyapunov_exponent();     
    }
     else if (eq == "Relay2damp") {
        // x' + b x = Sign( x_tau )
        double c = (double)json_params["c"];
        double gamma= (double)json_params["gamma"];
        double b0= (double)json_params["b0"];
        double b1= (double)json_params["b1"];
        double tau=(double)json_params["tau"];
        double t_finish = (double)json_params["t_finish"];
        double h = (double)json_params["h"];
        
        auto dde = Relay2damp(gamma, b0, b1, c, tau);
        Func<1,2> phi            = [&](double t) {return Vec<2>{0.5*d/c, 0};};
        Func<1,2> phi_derivative = [&](double t) {return Vec<2>{0, 0};};
        InitialValueProblem<2,1,0> ivp (dde, phi, phi_derivative);
        // ivp.compute_solution(h, t_finish);
        Func<1,2> psi = [tau](double t){return Vec<2>{sin(1.1*2*M_PI/tau*t), 1.1*2*M_PI/tau*cos(2*M_PI/tau*t)};};
        ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
        output = ivp.get_solution_and_lyapunov_exponent();     
    }
    else  if (eq == "RandomWalk") {
        // x' + b x = Sign( x_tau )
        double b = (double)json_params["b"];
        double d0= (double)json_params["d0"];
        double d1= (double)json_params["d1"];
        double tau=(double)json_params["tau"];
        double x0 = (double)json_params["x0"];
        double t_finish = (double)json_params["t_finish"];
        double h = (double)json_params["h"];
        
        auto dde = RandomWalk(b, d0, d1, tau);
        Func<1,1> phi            = [x0](double t) {return x0;};
        Func<1,1> phi_derivative = [](double t) {return 0;};
        InitialValueProblem<1,1,0> ivp (dde, phi, phi_derivative);
        // ivp.compute_solution(h, t_finish);
        Func<1,1> psi = [](double t){return 1.;};
        ivp.compute_solution_and_lyapunov_exponent(h, t_finish, psi);
        
        output = ivp.get_solution_and_lyapunov_exponent();
    } 
    
    
    
    string filename = output_prefix + " " + params;
    if (filename.size() > 200)
        filename.erase(200, string::npos);

	save(output, "output_bin/" + filename + ".bin");
    
	// save(output, "solution.bin");
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;

    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}
