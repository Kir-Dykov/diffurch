
/*
Implement:

    default.cpp
        here the user can externally specify the name of return handler and initial condition source
        usage: simple value, solution, dense solution, poincare sequence, etc
    
    batch of initial conditions -> array of dense output solutions (try this on double pendulum)
    
    1d:
        bifurcation diagram of a poincare section
    
    2d:
        period of a solution depending on parameters
        
    rk_test_vs_analytic

    rk_test
        here make it not dependent on the analytic solution
*/



// VARIABLE STEP SIZE!!!






























// how it should be

// Make solution
/*
    In python: (DONE!!!)
        from dde_solver import cpp_run
        params = dict(eq = "ODE_harmonic_oscillator", k = 3.)
        output = cpp_run("solution", params,  recompile=recompile, recalculate=recalculate, flags = "-O3")
        t = output[0]
        x, dx = output[1].T
        plt.plot(t, x)
        # plt.plot(x, dx)

    In cpp, you tell what to do:
        int main(int argc, char* argv[]) {
            json json_params = json::parse(argv[1]);
            string output_filename = argv[2];
            
            auto de            = de_from_json(json_params); // ???????????????????????????????????????????? from parameter "eq" the type of de is determied (at runtime)
            auto [t_finish, h] = unpack_json<double, "t_finish", "h">(json_params);
            auto phi = de.phi_parametrization_json(json_params); // ????
            save_arrays(output_file, de.solution<Return::All>(h, t_finish, phi));
            return 0;
        }
        
        
        int main(int argc, char* argv[]) {
            json json_params = json::parse(argv[1]);
            string output_filename = argv[2];
            
            auto de_constructor= make_de_constructor(json_params["eq"]);
            auto de            = de_constructor(json_params); // ????????????????????????????????????????????
            auto [t_finish, h] = unpack_json<double, "t_finish", "h">(json_params);
            auto phi = de.phi_parametrization_json(json_params); // ????
            save_arrays(output_file, de.solution<Return::All>(h, t_finish, phi));
            return 0;
        }
        
    In hpp, you define an equations:
        DDE_lin_1
            parameters = A, B
            variables = x
            f = x, x_tau -> A*x + B*x_tau)
            tau = 1
        
        ODE_lin_1
            parameters = k
            variables = x
            f = x -> k*x
            analytic_solution = t -> exp(k*t)
        
*/






// Make an image
/*
    In python: (DONE!!!)
        from dde_solver import cpp_run
        params = dict(x0 = 0)
        compiler_params = dict(eq = "ODE_harmonic_oscillator", X = k, Y = dx0)
        output = cpp_run("frequency", params, compiler_params, recompile=recompile, recalculate=recalculate, flags = "-O3")
        image = output.T

    In cpp, you tell what to do:
    
        // return last state after a time T
        int main(int argc, char* argv[]) {
            json json_params = json::parse(argv[1]);
            string output_filename = argv[2];
            
            JSON_UNPACK(json_params, double, X_l, X_r, Y_l, Y_r, t_finish, h);
            JSON_UNPACK(json_params, int, X_n, Y_n);
            
            auto X_s = linspace(X_l, X_r, X_n);
            auto Y_s = linspace(Y_l, Y_r, Y_n);
            array<array<double, Y_n>, X_n> image;
            
            auto de = EQ(json_params);
            auto phi = de.phi(json_params); // this can be different
            
            for (int i = 0; i < X_n; i++) {
                for (int j = 0; j < Y_n; j++) {
                    de.X = X_s[i];
                    de.Y = Y_s[j];
                    image[i][j] = de.solution<Return::Last>(h, t_finish, phi);
                }
            }
            
            save_arrays(output_filename, image);
            return 0;
        }
        
       
        
*/















//Ideas
/*
        UNPACK_PARAMETERS(double, json_params, alpha, tau); // macro that really does `auto [alpha, tau] = unpack_json<double>(json_params, "alpha", "tau");`

*/




/*
    Ranges of parameters
    
    tau --- single parameter
    tau_l tau_r tau_n


*/

/*
ode1(a, b, c)

for (i, j):
    a = as[i]
    b = bs[j]
    img[i,j] = F(ode(a,b, c))

*/

/*


auto params = json_unpack<EQ::params>(json_values)
auto de = EQ(params); 
// auto de = apply(EQ, params);




... -DVAR1=a -DVAR2=b

VAR1#_s = linspace(VAR1#_l,VAR1#_r,VAR1#_n);
VAR2#_s = linspace(VAR2#_l,VAR2#_r,VAR2#_n);

for (i, j):
    VAR1 = VAR1#_s[i]
    VAR2 = VAR2#_s[j]
    img[i,j] = F(ode(VAR1, VAR2, c))
    

*/


// COMMON USES
/*
    1. One set of parameters
        plot the solution, lyapunov exponents time series
    2. Iterating over some ranges of paramters
        D-section: graphing some property of a solution depending on various initial conditions or parameters

    1.
        auto params = json_unpack<EQ::params>(json_values)
        auto de = EQ(params); 
        // auto de = apply(EQ, params);
        
    2. Here I want to vary what parameters are to be varied, so the rest of parameters
        VAR1#_s = linspace(VAR1#_l,VAR1#_r,VAR1#_n);
        VAR2#_s = linspace(VAR2#_l,VAR2#_r,VAR2#_n);
        
        
        for (i, j):
            VAR1 = VAR1#_s[i]
            VAR2 = VAR2#_s[j]
            std::get<VAR1_index>(params) = VAR1#_s[i];
            std::get<VAR2_index>(params) = VAR2#_s[i];
            
            img[i,j] = F(ode(VAR1, VAR2, c))
*/