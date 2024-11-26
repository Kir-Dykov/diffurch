#include "../library/include.hpp"
using namespace std;

struct RelayNDDE1 : public DifferentialEquation<1, ArgSpec<true,  0, 0, 0>,  ArgSpec<false, 0, 1, 0>, 1> {
    // x' = k x + alpha sign(x_tau), alpha < 0
    
    // TODO make this into signle macro
    Real alpha, T;
    using param_names = ParamNames<"alpha", "T">;
    AUTO(param_tuple, tie(alpha, T));
    
    RelayNDDE1(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);   

        f[0] = {[this](const Vec<1>& x){
            return - x - alpha;
        }};
        f[1] = {[this](const Vec<1>& x){
            return - x + alpha;
        }};
        
        b = {[](Vec<1> dot_x_T){return dot_x_T;}};
        b_levels = {0.};
        
        b_spec = {{}, {Delay(T)}};  
    }
    
    auto analytic_periodic_solution(vector<Real> ts) {
        size_t k = ts.size() - 1;
        
        Real A = 0;
        
        Real pm = 1.;
        for (size_t j = 1; j <= k; j++) {
            pm *= -1.;
            A += pm*(exp(ts[j]) - exp(ts[j-1]));
        }
        A *= alpha * exp(-T);
        
        Real s = 1.;
        Real h_star = s*A/(minus_one_to_the(k) - exp(-T));
        vector<Real> Cs;
        Real c = (h_star + s*alpha);
        
        for (size_t i = 1; i <= k; i++) {
            Cs.push_back(c);
            c += 2*s*alpha*minus_one_to_the(i)*exp(ts[i]);
        }
        
        return VecMapC<1,1,1>{[this, ts, Cs, s](Real t) {
            t = fmod_positive(t, 2*T);
            Real s2 = 1.;
            if (t >= T) {
                t -= T;
                if (alpha < 0) {
                    s2 = -1.;
                }
            }
            
            int i = 1;
            while (ts[i] < t) i++;
            Real x = s*alpha*minus_one_to_the(i) + Cs[i-1]*exp(-t);
            return x*s2;
        }, {[this, ts, Cs, s](Real t) {
            t = fmod_positive(t, 2*T);
            Real s2 = 1.;
            if (t >= T) {
                t -= T;
                if (alpha < 0) {
                    s2 = -1.;
                }
            }
            
            int i = 1;
            while (ts[i] < t) i++;
            Real dot_x = -Cs[i-1]*exp(-t);
            return dot_x*s2;
        }}};
    }
};



struct RelayNDDE1arctan : public DifferentialEquation<1, ArgSpec<true,  0, 1, 0>,  ArgSpec<false, 0, 0, 0>, 0> {
    // x' = k x + alpha sign(x_tau), alpha < 0
    
    // TODO make this into signle macro
    Real alpha, T, eps;
    using param_names = ParamNames<"alpha", "T", "eps">;
    AUTO(param_tuple, tie(alpha, T, eps));
    
    RelayNDDE1arctan(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);   

        f[0] = {[this](const Vec<2>& X){
            auto& [x, dot_x_T] = X;
            return - x + alpha * 2.r / M_PI * atan(dot_x_T / eps);
        }};
        f[1] = {f[0]};
        
        
        f_spec = {{}, {Delay(T)}};  
        
//         b = {[](Vec<1> dot_x_tau){return dot_x_tau;}};
//         b_levels = {0.};
        
//         b_spec = {{}, {Delay(tau)}}; 
    }    
};




struct RelayNDDEexpsin : public DifferentialEquation<1, ArgSpec<true,  0, 1, 0>,  ArgSpec<false, 0, 0, 0>, 0> {
    Real alpha, T, eps, A, mu;
    int p;
    using param_names = ParamNames<"alpha", "T", "eps", "A", "mu", "p">;
    AUTO(param_tuple, tie(alpha, T, eps, A, mu, p));
    
    inline Real F(Real z) {
        return (exp(5*z)-1)/(exp(5*z)+1) + (A*sin(mu*pow(z,p)))/(z);
    }
    
    inline Real F_eps(Real z) {
        return F(z / eps);
    }
    
    RelayNDDEexpsin(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);   

        f[0] = {[this](const Vec<2>& X){
            auto& [x, dot_x_T] = X;
            return - x + alpha * F_eps(dot_x_T);
        }};
        f[1] = {f[0]};
        
        
        f_spec = {{}, {Delay(T)}};  
        
//         b = {[](Vec<1> dot_x_tau){return dot_x_tau;}};
//         b_levels = {0.};
        
//         b_spec = {{}, {Delay(tau)}}; 
    }    
};


struct RelayNDDEr1 : public DifferentialEquation<1, ArgSpec<true,  0, 1, 0>,  ArgSpec<false, 0, 0, 0>, 0> {
    Real alpha, T, eps, A, mu;
    int p;
    using param_names = ParamNames<"alpha", "T", "eps">;
    AUTO(param_tuple, tie(alpha, T, eps));
    
    inline Real F(Real z) {
        return (z+z*abs(z))/(1+z*z);
    }
    
    inline Real F_eps(Real z) {
        return F(z / eps);
    }
    
    RelayNDDEr1(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);   

        f[0] = {[this](const Vec<2>& X){
            auto& [x, dot_x_T] = X;
            return - x + alpha * F_eps(dot_x_T);
        }};
        f[1] = {f[0]};
        
        
        f_spec = {{}, {Delay(T)}};  
        
//         b = {[](Vec<1> dot_x_tau){return dot_x_tau;}};
//         b_levels = {0.};
        
//         b_spec = {{}, {Delay(tau)}}; 
    }    
};

struct RelayNDDE1arctan_sin : public DifferentialEquation<1, ArgSpec<true,  0, 1, 0>,  ArgSpec<false, 0, 0, 0>, 0> {
    Real alpha, T, eps, A, mu;
    int p;
    using param_names = ParamNames<"alpha", "T", "eps", "A", "mu", "p">;
    AUTO(param_tuple, tie(alpha, T, eps, A, mu, p));
    
    inline Real F(Real z) {
        return 2.r / M_PI * atan(z) + (z<=0 ? 0 : A * pow(sin(mu * z*z), p) / z);
    }
    
    inline Real F_eps(Real z) {
        return F(z / eps);
    }
    
    RelayNDDE1arctan_sin(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);   

        f[0] = {[this](const Vec<2>& X){
            auto& [x, dot_x_T] = X;
            return - x + alpha * F_eps(dot_x_T);
        }};
        f[1] = {f[0]};
        
        
        f_spec = {{}, {Delay(T)}};  
        
//         b = {[](Vec<1> dot_x_tau){return dot_x_tau;}};
//         b_levels = {0.};
        
//         b_spec = {{}, {Delay(tau)}}; 
    }    
};


struct RelayNDDEspikes : public DifferentialEquation<1, ArgSpec<true,  0, 1, 0>,  ArgSpec<false, 0, 0, 0>, 0> {
    Real alpha, T, eps, A, mu;
    int p;
    using param_names = ParamNames<"alpha", "T", "eps", "A", "mu", "p">;
    AUTO(param_tuple, tie(alpha, T, eps, A, mu, p));
    
    inline Real F(Real z) {
        return z <= 0 ? -1r + 1r/(1r - z) : 1r + (-1r + A*pow(sin(mu*z*z),p))/(1r + z);
    }
    
    inline Real F_eps(Real z) {
        return F(z / eps);
    }
    
    RelayNDDEspikes(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);   

        f[0] = {[this](const Vec<2>& X){
            auto& [x, dot_x_T] = X;
            return - x + alpha * F_eps(dot_x_T);
        }};
        f[1] = {f[0]};
        
        
        f_spec = {{}, {Delay(T)}};  
        
//         b = {[](Vec<1> dot_x_tau){return dot_x_tau;}};
//         b_levels = {0.};
        
//         b_spec = {{}, {Delay(tau)}}; 
    }    
};