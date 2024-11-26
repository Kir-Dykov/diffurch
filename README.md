# oddesa (Ordinary and Delay Differential Equations Solver and Analizer)

This is development project that aims to become python package and c++ library that implements efficient numerical methods for a variety of classes of differential equations, including:
* ordinary differential equations
* delay differential equations, including time- and state-dependent delays, as well as neutral type delays
* partial differential equations, that can be reduced to systems of ordinary and delay differential equations
* discontinuous systems of the above types

For such systems we compute:
* solutions, dense solutions
* poincare maps
* lyapunov exponents (including delayed and/or discontinuous systems!)
* periodicity of solutions

Also we provide easy ways that help to verify the actual order of convergence of numerical methods. 
This will allow the user to trust the numerical results.

# Why new numeric library?
As it usually happens, the motivation to write a new library lies in frustration with the existing tools:
* It is hard to verify the quality of numerical solutions.
  I want to have dedicated functions, that report the order
  of convergence of numerical methods in each case.
  This is the only way that ensures that the numerical method
  is appropriate for a problem and there are no bugs in implementation.
* I often study discontinuous sytems, so the numerical method must take that into account. Usually it is impossible or hard to control.
* I usually want to study equations for ranges of parameters, computiing bifurcation diagrams etc.
  It leads to massive computations, that require best performance and parallelization.
* I develop numerical methods for computation of Lyapunov exponents for discontinuous delayed systems.
  These developments are new, so I need an implementation of them. Also I don't know implementations
  for Lyapunov exponents of continuous delayed systems, since it is also relatively young area of study.
* I implemented numerical methods for concrete systems many times. And each time my code got increasingly complex,
  making it unmaintainable, forcing me to completely rewrite everything for every new equation I study.
  I want to hide that complexity behind some levels of abstraction, leaving per-case optimizations to the library and cpp compiler.

# Intended use

Pseudo code of intended usage:
```python
import numpy as np
import oddesa

lorenz_eq = oddesa.equation(
    # definition of the equation, that will be parsed into c++ code and just-in-time compiled
    """
        x' = sigma*(y - x); 
        y' = rho*x - y - x*z;
        z' = -beta*z + x*z;
    """,
    # definition of the variational eqaution, that is used to compute lyapunov exponents
    variational_equation="""
        _x' = sigma*(_y - _x);
        _y' = rho*x_ - y_ - x*z_ + x_*z;
        _z' = -beta*z_ + x_*z + x*z_;
    """
    # definitions of parameters, that appear in equation
    sigma = 10,
    beta = 8/3,
    rho=np.linspace(0, 28, 200), # if parameter is specified by an array, then solution will be computed for each value of parameter
    )

ls1 = lorenz.solution(interval=[0,10], ic=dict(x=0,y=0.1,z=0.2)) # solutions for each rho with the initial condition (x = 0, y = 0.1, z = 0.2)
ls2 = lorenz.solution(at=10, interval=[0,10], ic=np.linspace([0,0,0], [1,1,1], 100)) # solutions for each rho and each initial condition xyz-triple given in array linspace([0,0,0], [1,1,1], 100)
lorenz_le_s = lorenz.lyapunov_exponents(n="all", at="last", interval=[0, 10], ic=(0,0.1,0.2), variational_ic="random") // all 3 lyapunov exponents for each value of rho
```
