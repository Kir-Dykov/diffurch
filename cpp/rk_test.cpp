/*
I want to test my program for an order of a method.

The following systems seem appropriate:

ODE:
x' = x         x = e^t               x(0) = 1
x'' = -x       x = cos(t)           x(0) = 1, x'(0) = 0
x'' = x(1 - ln x - ln^2 x)         x = e^sin(t), x(0) = 1, x'(0) = 1

Discontinuous ODE:
x'' = -2sign(x)
x' = x + H(x - 1) * (x - 1)^q                 q is the order of discontinuity: q = 0 is just a jump, q = 1 is a jump of the derivative, etc

Continuous DDE
x' = x/2 + e/2 x(t - 1)     has solution x = C e^t
x' = e x(t - 1)             has solution x = C e^t
x' = - x ln(x(t-pi/2))      has solution x = e^{C sin(t)}

Discontinuous DDE
x' = -x + alpha sign(x_tau)     solution can be found explicitly

Continuous NDDE
x' = x'_tau // solution is known for periodic initial function

Discontinuous NDDE
x' + x = sign(x'_tau)

*/