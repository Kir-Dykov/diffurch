#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>

#include "library/json_unpack.hpp"
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




int main(int argc, char* argv[]) {

    // feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);


    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	string params = argv[1];
	json json_params = json::parse(params);
	cout << "~~~  parameters: " << params << " ~~~" << endl;
    
	string output_filename = argv[2];
    
    vector<vector<double>> output;
    
    auto [alpha, tau, t_finish, h, A, w] =
    unpack_json_doubles(json_params, array{"alpha", "tau", "t_finish", "h", "A", "w"});
        
    auto de = ClipDDE1(alpha, tau);
        
    VecFuncC<1,1> phi ([&](double t) {return Vec<1>{A*cos(w*t)};}, {[&](double t) {return Vec<1>{-w*A*sin(w*t)};}});
    
    output = de.solution<ReturnOption::AllInit>(h, t_finish, phi);
    
	save(output, "../output/bin/" + output_filename + ".bin");
    
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;

    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}