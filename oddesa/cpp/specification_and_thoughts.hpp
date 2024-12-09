#pragma once

#include "equation.hpp"

struct MyEQ : IVP<2, MyEq> {
    //parameter
    double k = 1;
    
    vector<function<double(double)>> delays {[this](double t){reutrn 1.}};
    double delay_max = 1.;
    
    function<Vec<n>(double, Vec<n>)> LHS = [this](double t, Vec<n> XXX) {
        auto& [x, y] = XXX;
        auto x|1  = state.evaluate<0>(t - delays[0](t))[0];
        auto x`|1 = state.evaluate<1>(t - delays[0](t))[1];
        return Vec<n>{
            -y - x|1 + k*x`|1,
            x|1
        };
    };
    
    

    
    function<double(Vec<n>)> initial_condition = [this](double t) {
        return Vec<n>{ sin(t), k*cos(t) };
    };
    
    function<double(Vec<n>)> initial_condition_derivative = [this](double t) {
        return Vec<n>{ cos(t), -k*sin(t) };
    };
    
    double integration_interval_start = 0;
    double integration_interval_finsih = 1;
    
    return_handler; ???
    
    stepsize_control; ???
        
    discrete_variables???
        
    tuple<vector<tuple<double>>> save;
        
    vector<Event> events = {
        Event {.detect = [&](){return X[0]*prevX[0] - 1 <= 0 && X[0] - 1 != 0;},
               .locate = [&](){return X[0] - 1;},
               .change = [&](){X[0] = -X[0];},
               .save   = [&](){save.get<0>().push_back(make_tuple(t));}
    };     
    
    auto solution() {...}
}

// can be fixed, or adaptive
// if fixed, all things here are trivial
// if adaptive, it requires error estimation and last stepsizes in State, order of RK
struct StepsizeController {
    double next_stepsize(State&);
    
    bool reject_last(State&);
    
    double corrected_last_stepsize(State&);
};

struct Delay {
    
}

struct Event {
    
}

struct State {
    
}

