#pragma once

#include "ddes.hpp"
#include "array_operations.hpp"




    

// x' + b x = Sign( x_tau )
RelayDDE<1, 1, 0> Relay1 (double b, double d0, double d1, double tau) {  
    const int n = 1;
    const int n_tau = 1;
    const int m_tau = 0;
    static const int n_arg_f = n*(1 + m_tau);
    using Func_b = Func<n*n_tau, 1>;
    using Func_db = Func<n*n_tau, n*n_tau>;
    using Func_f = Func<n_arg_f, n>;
    using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
    Func_b B = [](const double& X) {return X;};
    Func_db B_derivative = [](const double& X) {return 1.;};
    Func_f f0   =  [b,d0](const double& X) {return  -b*X + d0;};
    Func_f f1   =  [b,d1](const double& X) {return  -b*X + d1;};
    Func_df df0 =  [b,d0](const double& X, const double& W) {return -b*W;};
    Func_df df1 =  [b,d0](const double& X, const double& W) {return -b*W;};
    
    return RelayDDE<n, n_tau, m_tau>(B, B_derivative, f0, f1, df0, df1, tau);
}

// // x'' + b x' + c x = d sign( x_tau )
RelayDDE<2, 1, 0> Relay2 (double b, double c, double d, double tau) {
    const int n = 2;
    const int n_tau = 1;
    const int m_tau = 0;
    static const int n_arg_f = n*(1 + m_tau);
    using Func_b = Func<n*n_tau, 1>;
    using Func_db = Func<n*n_tau, n*n_tau>;
    using Func_f = Func<n_arg_f, n>;
    using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
    Func_b B = [](const Vec<n*n_tau>& X) {return X[0];};
    Func_db B_derivative = [](const Vec<n*n_tau>& X) {return Vec<2>{1., 0};};
    Func_f f0   = [b,c,d](const Vec<n_arg_f>& X) {return Vec<2>{X[1], -b*X[1] - c*X[0] - d};};
    Func_f f1   = [b,c,d](const Vec<n_arg_f>& X) {return Vec<2>{X[1], -b*X[1] - c*X[0] + d};};
    Func_df df0 = [b,c,d](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return Vec<2>{W[1], -b*W[1] - c*W[0]};};
    Func_df df1 = [b,c,d](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return Vec<2>{W[1], -b*W[1] - c*W[0]};};
    
    return RelayDDE<2, 1, 0>(B, B_derivative, f0, f1, df0, df1, tau);
}


// x' + b x_tau = Sign( x_tau )
RelayDDE<1, 1, 1> Croissant (double b, double d0, double d1, double tau) {  
    const int n = 1;
    const int n_tau = 1;
    const int m_tau = 1;
    static const int n_arg_f = n*(1 + m_tau);
    using Func_b = Func<n*n_tau, 1>;
    using Func_db = Func<n*n_tau, n*n_tau>;
    using Func_f = Func<n_arg_f, n>;
    using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
    Func_b B = [](const double& X) {return X;};
    Func_db B_derivative = [](const double& X) {return 1.;};
    Func_f f0   =  [b,d0](const Vec<n_arg_f>& X) {return  -b*X[1] + d0;};
    Func_f f1   =  [b,d1](const Vec<n_arg_f>& X) {return  -b*X[1] + d1;};
    Func_df df0 =  [b,d0](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return -b*W[1];};
    Func_df df1 =  [b,d0](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return -b*W[1];};
    
    return RelayDDE<n, n_tau, m_tau>(B, B_derivative, f0, f1, df0, df1, tau);
}

// x' + b x_tau = Sign( x_tau )
RelayDDE<1, 1, 1> MakeyGlass (double b, double d0, double d1, double tau) {  
    const int n = 1;
    const int n_tau = 1;
    const int m_tau = 1;
    static const int n_arg_f = n*(1 + m_tau);
    using Func_b = Func<n*n_tau, 1>;
    using Func_db = Func<n*n_tau, n*n_tau>;
    using Func_f = Func<n_arg_f, n>;
    using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
    Func_b B = [](const double& X) {return abs(X) - 1;};
    Func_db B_derivative = [](const double& X) {return sign(X);};
    Func_f f0   =  [b,d0](const Vec<n_arg_f>& X) {return  -b*X[0] + d0*X[1];};
    Func_f f1   =  [b,d1](const Vec<n_arg_f>& X) {return  -b*X[0] + d1*X[1];};
    Func_df df0 =  [b,d0,d1](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return -b*W[0] + d0*W[1];};
    Func_df df1 =  [b,d0,d1](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return -b*W[0] + d1*W[1];};
    
    return RelayDDE<n, n_tau, m_tau>(B, B_derivative, f0, f1, df0, df1, tau);
}

// x' + b x_tau = Sign( x_tau )
RelayDDE<1, 1, 1> MakeyGlassExp (double b, double d0, double d1, double tau) {  
    const int n = 1;
    const int n_tau = 1;
    const int m_tau = 1;
    static const int n_arg_f = n*(1 + m_tau);
    using Func_b = Func<n*n_tau, 1>;
    using Func_db = Func<n*n_tau, n*n_tau>;
    using Func_f = Func<n_arg_f, n>;
    using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
    Func_b B = [](const double& X) {return X;};
    Func_db B_derivative = [](const double& X) {return 1;};
    Func_f f0   =  [b,d0](const Vec<n_arg_f>& X) {return  -b + d0*exp(X[1] - X[0]);};
    Func_f f1   =  [b,d1](const Vec<n_arg_f>& X) {return  -b + d1*exp(X[1] - X[0]);};
    Func_df df0 =  [b,d0,d1](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return d0*exp(X[1] - X[0])*(W[1]-W[0]);};
    Func_df df1 =  [b,d0,d1](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return d1*exp(X[1] - X[0])*(W[1]-W[0]);};
    
    return RelayDDE<n, n_tau, m_tau>(B, B_derivative, f0, f1, df0, df1, tau);
}

// // x'' + b x' + c x = d sign( x_tau )
RelayDDE<2, 1, 0> Relay2damp (double gamma, double b0, double b1,  double c, double tau) {
    const int n = 2;
    const int n_tau = 1;
    const int m_tau = 0;
    static const int n_arg_f = n*(1 + m_tau);
    using Func_b = Func<n*n_tau, 1>;
    using Func_db = Func<n*n_tau, n*n_tau>;
    using Func_f = Func<n_arg_f, n>;
    using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
    Func_b B = [gamma](const Vec<n*n_tau>& X) {return X[1] - gamma;};
    Func_db B_derivative = [](const Vec<n*n_tau>& X) {return Vec<2>{1., 0};};
    Func_f f0   = [b0,c](const Vec<n_arg_f>& X) {return Vec<2>{X[1], -b0*X[1] - c*X[0]};};
    Func_f f1   = [b1,c](const Vec<n_arg_f>& X) {return Vec<2>{X[1], -b1*X[1] - c*X[0]};};
    Func_df df0 = [b0,c](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return Vec<2>{W[1], -b0*W[1] - c*W[0]};};
    Func_df df1 = [b1,c](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return Vec<2>{W[1], -b1*W[1] - c*W[0]};};
    
    return RelayDDE<2, 1, 0>(B, B_derivative, f0, f1, df0, df1, tau);
}

// x' + b x_tau = Sign( x_tau )
RelayDDE<1, 1, 0> RandomWalk (double b, double d0, double d1, double tau) {  
    const int n = 1;
    const int n_tau = 1;
    const int m_tau = 0;
    static const int n_arg_f = n*(1 + m_tau);
    using Func_b = Func<n*n_tau, 1>;
    using Func_db = Func<n*n_tau, n*n_tau>;
    using Func_f = Func<n_arg_f, n>;
    using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
    Func_b B = [](const double& X) {return sin(X);};
    Func_db B_derivative = [](const double& X) {return cos(X);};
    Func_f f0   =  [b,d0](const Vec<n_arg_f>& X) {return  -b*X + d0;};
    Func_f f1   =  [b,d1](const Vec<n_arg_f>& X) {return  -b*X + d1;};
    Func_df df0 =  [b,d0](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return -b*W;};
    Func_df df1 =  [b,d0](const Vec<n_arg_f>& X, const Vec<n_arg_f>& W) {return -b*W;};
    
    return RelayDDE<n, n_tau, m_tau>(B, B_derivative, f0, f1, df0, df1, tau);
}
    

// // x' + eps x = d Sign( sin(x_tau) )
// RelayDDE<1, 1, 0> RandomWalk (double eps, double d0, double d1, double tau) {
//     const int n = 1;
//     auto B = [](const array<double, n>& X) {return sin(X[0]);};
//     auto B_derivative = [](const array<double, n>& X) {return array<double, n>{cos(X[0])};};
//     auto f0 =[eps,d0](const array<double, n>& X) {return array<double, n>{ -eps*X[0] + d0};};
//     auto f1 =[eps,d1](const array<double, n>& X) {return array<double, n>{ -eps*X[0] + d1};};
    
//     return RelayDDE<n, 1, 0>(B, B_derivative, f0, f1, {tau});
// }

// // x' + b x_tau = d Sign( x_tau )
// RelayDDE<1, 1, 1> Croissant (double b, double d0, double d1, double tau) {
//     const int n = 1;
//     auto B = [](const array<double, n>& X) {return X[0];};
//     auto B_derivative = [](const array<double, n>& X) {return array<double, n>{1};};
//     auto f0 =[b,d0](const array<double, n*2>& X) {return array<double, n>{ -b*X[1] + d0};};
//     auto f1 =[b,d1](const array<double, n*2>& X) {return array<double, n>{ -b*X[1] + d1};};
//     return RelayDDE<n, 1, 1>(B, B_derivative, f0, f1, {tau});
// }

// // x' + b x = d x_tau H(1-|x_tau|)
// RelayDDE<1, 1, 1> LimitMakeyGlass (double b, double d, double tau) {
//     const int n = 1;
//     auto B = [](const array<double, n>& X) {return abs(X[0]) - 1;};
//     auto B_derivative = [](const array<double, n>& X) {return array<double, n>{ sign(X[0] )};};
//     auto f0 =[b,d](const array<double, n*2>& X) {return array<double, n>{ -b*X[0] + d*X[1]};};
//     auto f1 =[b](const array<double, n*2>& X) {return array<double, n>{ -b*X[0]};};
//     return RelayDDE<n, 1, 1>(B, B_derivative, f0, f1, {tau});
// }

// // 
// RelayDDE<2, 1, 0> OpenBeak (double nu, double gamma, double d0, double d1, double tau) {
//     const int n = 2;
//     auto B = [gamma](const array<double, n>& X) {return gamma*gamma*X[0]*X[0] + X[1]*X[1] - 1;};
//     auto B_derivative = [b,c,d](const array<double, n>& X) {return array<double, n>{2*gamma*gamma*X[0], 2*X[1]};};
//     auto f0 =[nu,d0,d1](const array<double, n>& X) {return array<double, n>{nu*X[1] + d0*X[0], -nu*X[0] + d0*X[1]};};
//     auto f1 =[nu,d0,d1](const array<double, n>& X) {return array<double, n>{nu*X[1] + d1*X[0], -nu*X[0] + d1*X[1]};};
//     return RelayDDE<n, 1, 0>(B, B_derivative, f0, f1, {tau});
// }



// // 
// RelayDDE<2, 1, 0> DampingControl (double gamma, double d0, double d1, double tau) {
//     const int n = 2;
    
//     auto B = [b,c,d](const array<double, n>& X) {return X[1] - gamma;};
//     auto B_derivative = [b,c,d](const array<double, n>& X) {return array<double, n>{1, 0};};
//     auto f0 =[b,c,d](const array<double, n>& X) {return array<double, n>{X[1], -d0*X[1] - X[0]};};
//     auto f1 =[b,c,d](const array<double, n>& X) {return array<double, n>{X[1], -d1*X[1] - X[0]};};
//     return RelayDDE<n, 1, 0>(B, B_derivative, f0, f1, {tau});
// }