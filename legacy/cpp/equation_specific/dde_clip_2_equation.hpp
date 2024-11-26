#include "../library/include.hpp"
using namespace std;

struct DDE_clip_2 : public DifferentialEquation<2, ArgSpec<true,  1, 0, 0>,  ArgSpec<false, 1, 0, 0>, 2> {
    // x' = k x + alpha sign(x_tau), alpha < 0
    double alpha, sigma, tau;
    using param_names = ParamNames<"alpha", "sigma", "tau">;
    AUTO(param_tuple, tie(alpha, sigma, tau));
    
    DDE_clip_2(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);   

        f[0] = {[this](const Vec<4>& X){
            auto& [x, dot_x, x_tau, dot_x_tau] = X;
            return Vec<2>{dot_x, -sigma*dot_x - x - alpha};
        }};
        f[1] = {[this](const Vec<4>& X){
            auto& [x, dot_x, x_tau, dot_x_tau] = X;
            return Vec<2>{dot_x, -sigma*dot_x - x + alpha * x_tau};
        }};
        f[2] =  {[this](const Vec<4>& X){
            auto& [x, dot_x, x_tau, dot_x_tau] = X;
            return Vec<2>{dot_x, -sigma*dot_x - x + alpha};
        }};
        
        b = {[](Vec<2> X_tau){return X_tau[0];}};
        b_levels = {-1., 1.};
        
        f_spec = {{Delay(tau)}, {}};
        b_spec = {{Delay(tau)}, {}};
    

    }
};