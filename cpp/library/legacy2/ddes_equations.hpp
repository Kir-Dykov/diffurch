#pragma once

#include "dde.hpp"
#include "vec.hpp"

/* *******************************

This file contains the definitions of some concreate equations.
This is done in the form of the functions that take parameters 
that a particular equation contains and returns a DDE object
that corresponds to the equation with those parameters.

******************************* */

// x' + b x = d0 + (d1 - d0) H( x_tau )

// ArgSpec<true, 2, 2> // for f(x, x_tau1, x_tau2, x'_tau1, x'_tau2)
// ArgSpec<n, false, 3>// for b(x_tau1, x_tau2, x_tau3)
// b_levels

struct SignDDE1 : public DE<1, ArgSpec<true,  0, 0, 1>,  ArgSpec<false, 1, 0, 1>, 1> {
    SignDDE1(Real c, Real d0, Real d1, Real tau) {
        f = {
            VecMapC<1,1,1>([c,d0,d1](Real x){return -c*x + d0;}, {[c,d0,d1](Real x, Real delta_x) {return -c*delta_x;}}),
            VecMapC<1,1,1>([c,d0,d1](Real x){return -c*x + d1;}, {[c,d0,d1](Real x, Real delta_x) {return -c*delta_x;}}),
        };
        b = VecMapC<1,1,1>([](Real x){return x;}, {[](Real x, Real delta_x){return delta_x;}});
        
        b_levels = {0};
        b_delays = {tau};
        
        // f =  {[c,d0,d1](const Vec<1>& X) {return  Vec<1>{-c*X[0] + d0};},
        //       [c,d0,d1](const Vec<1>& X) {return  Vec<1>{-c*X[0] + d1};}};
        // df = {[c,d0,d1](const Vec<1>& X, const Vec<1>& W) {return Vec<1>{-c*W[0]};},
        //       [c,d0,d1](const Vec<1>& X, const Vec<1>& W) {return Vec<1>{-c*W[0]};}};
        // b = []      (const Vec<1>& X)                  {return X[0];};
        // b_derivative = [](const Vec<1>& X) {return Vec<1>{1.};};
        // b_levels = {0};
        // b_delays = {tau};
    }   
};

struct SignDDE2 : public DE<2, ArgSpec<true,  0, 0, 1>,  ArgSpec<false, 1, 0, 1>, 1> {
    SignDDE2(Real b_, Real c, Real d, Real tau) {
        f = {
            VecMapC<2,2,1>( [b_,c,d](const Vec<2>& X)                  {return  Vec<2>{X[1], -b_*X[1] - c*X[0] - d};}, 
                    {[b_,c,d](const Vec<2>& X, const Vec<2>& W) {return  Vec<2>{W[1], -b_*W[1] - c*W[0]};  }}),
            VecMapC<2,2,1>( [b_,c,d](const Vec<2>& X)                  {return  Vec<2>{X[1], -b_*X[1] - c*X[0] + d};}, 
                    {[b_,c,d](const Vec<2>& X, const Vec<2>& W) {return  Vec<2>{W[1], -b_*W[1] - c*W[0]};  }}),
        };
        
        b = VecMapC<2,1,1>([](const Vec<2>& X){return X[0];}, {[](const Vec<2>& X, const Vec<2>& W){return W[0];}});
        
        
        // f = { [b_,c,d](const Vec<2>& X) {return  Vec<2>{X[1], -b_*X[1] - c*X[0] - d};},
        //       [b_,c,d](const Vec<2>& X) {return  Vec<2>{X[1], -b_*X[1] - c*X[0] + d};}};
        // df = {[b_,c,d](const Vec<2>& X, const Vec<2>& W) {return  Vec<2>{W[1], -b_*W[1] - c*W[0]};},
        //       [b_,c,d](const Vec<2>& X, const Vec<2>& W) {return  Vec<2>{W[1], -b_*W[1] - c*W[0]};}};
        // b  = VecMacC([](const Vec<2>& X){return X[0];}, {[](const Vec<2>& X, const Vec<2>& W){return W[0];}});

        b_levels = {0};
        b_delays = {tau};
    }
};

struct ClipDDE1 : public  DE<1, ArgSpec<true,  1, 0, 1>,  ArgSpec<false, 1, 0, 1>, 2> {
    ClipDDE1(Real alpha, Real tau) {
        
        f = {
            VecMapC<2,1,1>([alpha](Vec<2> x){return -x[0] - alpha;},      {[alpha](Vec<2> x, Vec<2> delta_x) {return -delta_x[0];}}),
            VecMapC<2,1,1>([alpha](Vec<2> x){return -x[0] + alpha*x[1];}, {[alpha](Vec<2> x, Vec<2> delta_x) {return -delta_x[0] + alpha*delta_x[1];}}),
            VecMapC<2,1,1>([alpha](Vec<2> x){return -x[0] + alpha;},      {[alpha](Vec<2> x, Vec<2> delta_x) {return -delta_x[0];}})
        };
        b = VecMapC<1,1,1>([](Real x){return x;}, {[](Real x, Real delta_x){return delta_x;}});
        
        b_levels = {-1, 1};
        
        f_delays = {tau};
        b_delays = {tau};
    }   
};






// // // x'' + b x' + c x = d sign( x_tau )
// DDE<2, 1, 0> Relay2 (Real b, Real c, Real d, Real tau) {
//     const int n = 2;
//     const int n_tau = 1;
//     const int m_tau = 0;
//     static const int n_arg_f = n*(1 + m_tau);
//     using Func_b = Func<n*n_tau, 1>;
//     using Func_db = Func<n*n_tau, n*n_tau>;
//     using Func_f = Func<n_arg_f, n>;
//     using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
//     Func_b B = [](const Vec<n*n_tau>& X) {return X[0];};
//     Func_db B_derivative = [](const Vec<n*n_tau>& X) {return Vec<2>{1., 0};};
//     Func_f f0   = [b,c,d](const Vec<n_arg_f>& X) {return Vec<2>{X[1], -b*X[1] - c*X[0] - d};};
//     Func_f f1   = [b,c,d](const Vec<n_arg_f>& X) {return Vec<2>{X[1], -b*X[1] - c*X[0] + d};};
//     Func_df df0 = [b,c,d](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return Vec<2>{W[1], -b*W[1] - c*W[0]};};
//     Func_df df1 = [b,c,d](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return Vec<2>{W[1], -b*W[1] - c*W[0]};};
    
//     return DDE<2, 1, 0>(B, B_derivative, f0, f1, df0, df1, tau);
// }


// // x' + b x_tau = Sign( x_tau )
// DDE<1, 1, 1> Croissant (Real b, Real d0, Real d1, Real tau) {  
//     const int n = 1;
//     const int n_tau = 1;
//     const int m_tau = 1;
//     static const int n_arg_f = n*(1 + m_tau);
//     using Func_b = Func<n*n_tau, 1>;
//     using Func_db = Func<n*n_tau, n*n_tau>;
//     using Func_f = Func<n_arg_f, n>;
//     using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
//     Func_b B = [](const Real& X) {return X;};
//     Func_db B_derivative = [](const Real& X) {return 1.;};
//     Func_f f0   =  [b,d0](const Vec<n_arg_f>& X) {return  -b*X[1] + d0;};
//     Func_f f1   =  [b,d1](const Vec<n_arg_f>& X) {return  -b*X[1] + d1;};
//     Func_df df0 =  [b,d0](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return -b*W[1];};
//     Func_df df1 =  [b,d0](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return -b*W[1];};
    
//     return DDE<n, n_tau, m_tau>(B, B_derivative, f0, f1, df0, df1, tau);
// }

// // x' + b x_tau = Sign( x_tau )
// DDE<1, 1, 1> MakeyGlass (Real b, Real d0, Real d1, Real tau) {  
//     const int n = 1;
//     const int n_tau = 1;
//     const int m_tau = 1;
//     static const int n_arg_f = n*(1 + m_tau);
//     using Func_b = Func<n*n_tau, 1>;
//     using Func_db = Func<n*n_tau, n*n_tau>;
//     using Func_f = Func<n_arg_f, n>;
//     using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
//     Func_b B = [](const Real& X) {return abs(X) - 1;};
//     Func_db B_derivative = [](const Real& X) {return sign(X);};
//     Func_f f0   =  [b,d0](const Vec<n_arg_f>& X) {return  -b*X[0] + d0*X[1];};
//     Func_f f1   =  [b,d1](const Vec<n_arg_f>& X) {return  -b*X[0] + d1*X[1];};
//     Func_df df0 =  [b,d0,d1](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return -b*W[0] + d0*W[1];};
//     Func_df df1 =  [b,d0,d1](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return -b*W[0] + d1*W[1];};
    
//     return DDE<n, n_tau, m_tau>(B, B_derivative, f0, f1, df0, df1, tau);
// }

// // x' + b x_tau = Sign( x_tau )
// DDE<1, 1, 1> MakeyGlassExp (Real b, Real d0, Real d1, Real tau) {  
//     const int n = 1;
//     const int n_tau = 1;
//     const int m_tau = 1;
//     static const int n_arg_f = n*(1 + m_tau);
//     using Func_b = Func<n*n_tau, 1>;
//     using Func_db = Func<n*n_tau, n*n_tau>;
//     using Func_f = Func<n_arg_f, n>;
//     using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
//     Func_b B = [](const Real& X) {return X;};
//     Func_db B_derivative = [](const Real& X) {return 1;};
//     Func_f f0   =  [b,d0](const Vec<n_arg_f>& X) {return  -b + d0*exp(X[1] - X[0]);};
//     Func_f f1   =  [b,d1](const Vec<n_arg_f>& X) {return  -b + d1*exp(X[1] - X[0]);};
//     Func_df df0 =  [b,d0,d1](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return d0*exp(X[1] - X[0])*(W[1]-W[0]);};
//     Func_df df1 =  [b,d0,d1](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return d1*exp(X[1] - X[0])*(W[1]-W[0]);};
    
//     return DDE<n, n_tau, m_tau>(B, B_derivative, f0, f1, df0, df1, tau);
// }

// // // x'' + b x' + c x = d sign( x_tau )
// DDE<2, 1, 0> Relay2damp (Real gamma, Real b0, Real b1,  Real c, Real tau) {
//     const int n = 2;
//     const int n_tau = 1;
//     const int m_tau = 0;
//     static const int n_arg_f = n*(1 + m_tau);
//     using Func_b = Func<n*n_tau, 1>;
//     using Func_db = Func<n*n_tau, n*n_tau>;
//     using Func_f = Func<n_arg_f, n>;
//     using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
//     Func_b B = [gamma](const Vec<n*n_tau>& X) {return X[1] - gamma;};
//     Func_db B_derivative = [](const Vec<n*n_tau>& X) {return Vec<2>{1., 0};};
//     Func_f f0   = [b0,c](const Vec<n_arg_f>& X) {return Vec<2>{X[1], -b0*X[1] - c*X[0]};};
//     Func_f f1   = [b1,c](const Vec<n_arg_f>& X) {return Vec<2>{X[1], -b1*X[1] - c*X[0]};};
//     Func_df df0 = [b0,c](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return Vec<2>{W[1], -b0*W[1] - c*W[0]};};
//     Func_df df1 = [b1,c](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return Vec<2>{W[1], -b1*W[1] - c*W[0]};};
    
//     return DDE<2, 1, 0>(B, B_derivative, f0, f1, df0, df1, tau);
// }

// // x' + b x_tau = Sign( x_tau )
// DDE<1, 1, 0> RandomWalk (Real b, Real d0, Real d1, Real tau) {  
//     const int n = 1;
//     const int n_tau = 1;
//     const int m_tau = 0;
//     static const int n_arg_f = n*(1 + m_tau);
//     using Func_b = Func<n*n_tau, 1>;
//     using Func_db = Func<n*n_tau, n*n_tau>;
//     using Func_f = Func<n_arg_f, n>;
//     using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
//     Func_b B = [](const Real& X) {return sin(X);};
//     Func_db B_derivative = [](const Real& X) {return cos(X);};
//     Func_f f0   =  [b,d0](const Vec<n_arg_f>& X) {return  -b*X + d0;};
//     Func_f f1   =  [b,d1](const Vec<n_arg_f>& X) {return  -b*X + d1;};
//     Func_df df0 =  [b,d0](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return -b*W;};
//     Func_df df1 =  [b,d0](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return -b*W;};
    
//     return DDE<n, n_tau, m_tau>(B, B_derivative, f0, f1, df0, df1, tau);
// }
    

// // x' + eps x = d Sign( sin(x_tau) )
// DDE<1, 1, 0> RandomWalk (Real eps, Real d0, Real d1, Real tau) {
//     const int n = 1;
//     auto B = [](const array<Real, n>& X) {return sin(X[0]);};
//     auto B_derivative = [](const array<Real, n>& X) {return array<Real, n>{cos(X[0])};};
//     auto f0 =[eps,d0](const array<Real, n>& X) {return array<Real, n>{ -eps*X[0] + d0};};
//     auto f1 =[eps,d1](const array<Real, n>& X) {return array<Real, n>{ -eps*X[0] + d1};};
    
//     return DDE<n, 1, 0>(B, B_derivative, f0, f1, {tau});
// }

// // x' + b x_tau = d Sign( x_tau )
// DDE<1, 1, 1> Croissant (Real b, Real d0, Real d1, Real tau) {
//     const int n = 1;
//     auto B = [](const array<Real, n>& X) {return X[0];};
//     auto B_derivative = [](const array<Real, n>& X) {return array<Real, n>{1};};
//     auto f0 =[b,d0](const array<Real, n*2>& X) {return array<Real, n>{ -b*X[1] + d0};};
//     auto f1 =[b,d1](const array<Real, n*2>& X) {return array<Real, n>{ -b*X[1] + d1};};
//     return DDE<n, 1, 1>(B, B_derivative, f0, f1, {tau});
// }

// // x' + b x = d x_tau H(1-|x_tau|)
// DDE<1, 1, 1> LimitMakeyGlass (Real b, Real d, Real tau) {
//     const int n = 1;
//     auto B = [](const array<Real, n>& X) {return abs(X[0]) - 1;};
//     auto B_derivative = [](const array<Real, n>& X) {return array<Real, n>{ sign(X[0] )};};
//     auto f0 =[b,d](const array<Real, n*2>& X) {return array<Real, n>{ -b*X[0] + d*X[1]};};
//     auto f1 =[b](const array<Real, n*2>& X) {return array<Real, n>{ -b*X[0]};};
//     return DDE<n, 1, 1>(B, B_derivative, f0, f1, {tau});
// }

// // 
// DDE<2, 1, 0> OpenBeak (Real nu, Real gamma, Real d0, Real d1, Real tau) {
//     const int n = 2;
//     auto B = [gamma](const array<Real, n>& X) {return gamma*gamma*X[0]*X[0] + X[1]*X[1] - 1;};
//     auto B_derivative = [b,c,d](const array<Real, n>& X) {return array<Real, n>{2*gamma*gamma*X[0], 2*X[1]};};
//     auto f0 =[nu,d0,d1](const array<Real, n>& X) {return array<Real, n>{nu*X[1] + d0*X[0], -nu*X[0] + d0*X[1]};};
//     auto f1 =[nu,d0,d1](const array<Real, n>& X) {return array<Real, n>{nu*X[1] + d1*X[0], -nu*X[0] + d1*X[1]};};
//     return DDE<n, 1, 0>(B, B_derivative, f0, f1, {tau});
// }



// // 
// DDE<2, 1, 0> DampingControl (Real gamma, Real d0, Real d1, Real tau) {
//     const int n = 2;
    
//     auto B = [b,c,d](const array<Real, n>& X) {return X[1] - gamma;};
//     auto B_derivative = [b,c,d](const array<Real, n>& X) {return array<Real, n>{1, 0};};
//     auto f0 =[b,c,d](const array<Real, n>& X) {return array<Real, n>{X[1], -d0*X[1] - X[0]};};
//     auto f1 =[b,c,d](const array<Real, n>& X) {return array<Real, n>{X[1], -d1*X[1] - X[0]};};
//     return DDE<n, 1, 0>(B, B_derivative, f0, f1, {tau});
// }