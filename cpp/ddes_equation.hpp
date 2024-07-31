#pragma_once

#include "ddes.hpp"

// x'' + b x' + c x = d sign( x_tau )
RelayDDE<2, 1, 0> Relay2 (double b, double c, double d, double tau) {
    const int n = 2;
    
    auto B = [b,c,d](const array<double, n>& X) {return X[0];};
    auto B_derivative = [b,c,d](const array<double, n>& X) {return array<double, n>{1, 0};};
    auto f0 =[b,c,d](const array<double, n>& X) {return array<double, n>{X[1], -b*X[1] - c*X[0] - d};};
    auto f1 =[b,c,d](const array<double, n>& X) {return array<double, n>{X[1], -b*X[1] - c*X[0] + d};};
    return RelayDDE<n, 1, 0>(B, B_derivative, f0, f1, {tau});
}
    

// x' + b x = Sign( x_tau )
RelayDDE<1, 1, 0> Relay1 (double b, double d0, double d1, double tau) {
    const int n = 1;
    auto B = [](const array<double, n>& X) {return X[0];};
    auto B_derivative = [](const array<double, n>& X) {return array<double, n>{1};};
    auto f0 =[b,d0](const array<double, n>& X) {return array<double, n>{ -b*X[0] + d0};};
    auto f1 =[b,d1](const array<double, n>& X) {return array<double, n>{ -b*X[0] + d1};};
    
    auto df0 =  [b,d0](const array<double, n>& X, const array<double, n>& W) {return array<double, n>{ -b*W[0]};};
    auto df1 =  [b,d0](const array<double, n>& X, const array<double, n>& W) {return array<double, n>{ -b*W[0]};};
    
    return RelayDDE<n, 1, 0>(B, B_derivative, f0, f1, {tau});
}
    

// x' + eps x = d Sign( sin(x_tau) )
RelayDDE<1, 1, 0> RandomWalk (double eps, double d0, double d1, double tau) {
    const int n = 1;
    auto B = [](const array<double, n>& X) {return sin(X[0]);};
    auto B_derivative = [](const array<double, n>& X) {return array<double, n>{cos(X[0])};};
    auto f0 =[eps,d0](const array<double, n>& X) {return array<double, n>{ -eps*X[0] + d0};};
    auto f1 =[eps,d1](const array<double, n>& X) {return array<double, n>{ -eps*X[0] + d1};};
    
    return RelayDDE<n, 1, 0>(B, B_derivative, f0, f1, {tau});
}

// x' + b x_tau = d Sign( x_tau )
RelayDDE<1, 1, 1> Croissant (double b, double d0, double d1, double tau) {
    const int n = 1;
    auto B = [](const array<double, n>& X) {return X[0];};
    auto B_derivative = [](const array<double, n>& X) {return array<double, n>{1};};
    auto f0 =[b,d0](const array<double, n*2>& X) {return array<double, n>{ -b*X[1] + d0};};
    auto f1 =[b,d1](const array<double, n*2>& X) {return array<double, n>{ -b*X[1] + d1};};
    return RelayDDE<n, 1, 1>(B, B_derivative, f0, f1, {tau});
}

// x' + b x = d x_tau H(1-|x_tau|)
RelayDDE<1, 1, 1> LimitMakeyGlass (double b, double d, double tau) {
    const int n = 1;
    auto B = [](const array<double, n>& X) {return abs(X[0]) - 1;};
    auto B_derivative = [](const array<double, n>& X) {return array<double, n>{ sign(X[0] )};};
    auto f0 =[b,d](const array<double, n*2>& X) {return array<double, n>{ -b*X[0] + d*X[1]};};
    auto f1 =[b](const array<double, n*2>& X) {return array<double, n>{ -b*X[0]};};
    return RelayDDE<n, 1, 1>(B, B_derivative, f0, f1, {tau});
}

// 
RelayDDE<2, 1, 0> OpenBeak (double nu, double gamma, double d0, double d1, double tau) {
    const int n = 2;
    auto B = [gamma](const array<double, n>& X) {return gamma*gamma*X[0]*X[0] + X[1]*X[1] - 1;};
    auto B_derivative = [b,c,d](const array<double, n>& X) {return array<double, n>{2*gamma*gamma*X[0], 2*X[1]};};
    auto f0 =[nu,d0,d1](const array<double, n>& X) {return array<double, n>{nu*X[1] + d0*X[0], -nu*X[0] + d0*X[1]};};
    auto f1 =[nu,d0,d1](const array<double, n>& X) {return array<double, n>{nu*X[1] + d1*X[0], -nu*X[0] + d1*X[1]};};
    return RelayDDE<n, 1, 0>(B, B_derivative, f0, f1, {tau});
}



// 
RelayDDE<2, 1, 0> DampingControl (double gamma, double d0, double d1, double tau) {
    const int n = 2;
    
    auto B = [b,c,d](const array<double, n>& X) {return X[1] - gamma;};
    auto B_derivative = [b,c,d](const array<double, n>& X) {return array<double, n>{1, 0};};
    auto f0 =[b,c,d](const array<double, n>& X) {return array<double, n>{X[1], -d0*X[1] - X[0]};};
    auto f1 =[b,c,d](const array<double, n>& X) {return array<double, n>{X[1], -d1*X[1] - X[0]};};
    return RelayDDE<n, 1, 0>(B, B_derivative, f0, f1, {tau});
}