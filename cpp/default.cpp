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

#include "library/include.hpp"

using namespace std;

#ifndef RETURN
#define ReturnSolution()
#endif

#define MEASURE_TIME_BEGIN  std::chrono::steady_clock::time_point measure_time_begin = std::chrono::steady_clock::now();
#define MEASURE_TIME_END    std::chrono::steady_clock::time_point measure_time_end = chrono::steady_clock::now();
#define MEASURE_TIME_PRINT  int seconds = chrono::duration_cast<chrono::seconds>(measure_time_end - measure_time_begin).count(); \
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;


int main(int argc, char* argv[]) {
    FILE_IS_EXECUTED;
    MEASURE_TIME_BEGIN;
    
    json   JSON_PARAMS     = json::parse(argv[1]);
	string output_filename = argv[2];
    
    JSON_UNPACK(Real, t_finish, h);    
    auto de = from_json<EQ>(JSON_PARAMS);   
    auto ret = de.solution(h, t_finish, INIT_FUNCTION, RETURN);
    
    string full_output_filename = "../output/bin/" + output_filename + ".bin";
    
    if constexpr (is_tuple<decltype(ret)>::value) {
        std::apply([full_output_filename](auto&&... args){
            save_arrays(full_output_filename, std::forward<decltype(args)>(args)...);
        }, ret);
    } else {
        save_arrays(full_output_filename, ret);
    }
    
    MEASURE_TIME_END;
    MEASURE_TIME_PRINT;
    FILE_IS_FINISHED;
}



