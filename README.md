# oddesa (Ordinary and Delay Differential Equations Solver and Analyzer)

This is development project that aims to become python package and c++ library that implements efficient numerical methods for a variety of classes of differential equations, including:
* ordinary differential equations
* delay differential equations, including time- and state-dependent delays, as well as neutral type delays
* partial differential equations, that can be reduced to systems of ordinary and delay differential equations
* discontinuous systems of the above types

For such systems we compute:
* solutions, dense solutions
* Poincare maps
* Lyapunov exponents (including delayed and/or discontinuous systems!)
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
* I often study discontinuous systems, so the numerical method must take that into account. Usually it is impossible or hard to control.
* I usually want to study equations for ranges of parameters, computing bifurcation diagrams etc.
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

# Specification for future implementation
```python
oddesa.solve(equation, 
             events = None, 
             ic = None, 
             interval = None,
             at = "all",
             return_derivatives = 0,
             return_delayed = None,
             stepsize = ("auto", dict(tol =1e-10, atol = None, rtol = None),
             max_stepsize = 1.,
             max_steps = 0,
             events = None,
             disable_discontinuity_detection = True,
             discontinuities = None
)

```
Parameters: 
* `equation : str` 
  defines the system of differential equations by string in which:
  - equations are separated by `;` (semicolon)
  - each equation has the form `<variable>' = <rhs>`, 
    each equation must be solved for the highest order derivative
  - `<variable>` has the form of correct variable identifier 
    followed by any number of `'` (single quote), for example
    equation `x'' = x' + exp(-t)` by order reduction expands 
    into `(x)' = x'; (x')' = x' + exp(-t)` where the variables are `x'` and `x` (this one is implicit)
  - `<rhs>` defines the right hand side function 
    in terms of variables and parameters;
    all identifiers that were not recognized 
    as variables or known functions (like `exp`, `sin`, etc.)
    are interpreted as parameters
  - variables in `<rhs>` can appear by themselves or 
    in the form `<variable>|<delay>` denoting delayed variable `<variable>(t - <delay>)` 
    or `<variable>'|<delay>` denoting neutral delayed terms `<variable>'(t - <delay>)`
  - `<delay>` is either positive numeric literal or a single identifier, examples: `x|1`, `x|1.5`, `x|tau`.
  examples:
  - `x' = -k*x` : exponential decay
  - `x'' = - mu * x' + k*k * x` : dumped harmonic oscillator
  - `x'' = k_x*x + gamma*(y - x); y'' = k_y*y + gamma*(x - y)` : coupled harmonic oscillators
  - `x' = a*x + b*x|tau` : linear delay equation
  - `x' = a*x + b*x|tau + c*x'|tau` : linear neutral delay equation
* `ic : str | tuple of floats | tuple of arrays | array of arrays | dict of floats | dict of arrays` 
  defines the initial conditions
  
* `interval : (float, float)`
  defines the interval of integration
* `at : "none" | "all" | "last" | "finish" | None | array[float] | int`
  defines the values at which the solution is evaluated:
  - `"none" or None` : no integration points are returned.
    This option is useful when only output from event functions is needed
  - `"last" | "finish"` : return only the value at the end of integration
  - `"all"`  : the solution is evaluated at all integration points
  - `array[float]` : the solution is evaluated at the points specified in that array,
    points that are not in the `interval` are ignored, points after `stop_integration` event are ignored
  - `int` : the solution is evaluated at equally spaced points dividing the `interval`, 
    number of which is specified by this parameter, usage `at = int(N)` 
    is equivalent to `at = np.linspace(interval[0], interval[1], int(N))`, 
    the subdivision process is not affected by `stop_integration` event, 
    but the points beyond `stop_integration` are not evaluated
* `return_derivatives : int | None | array_like `
  the number of derivatives of solutions evaluated
  - `None | 0` : return only the solution, don't evaluate it's derivatives
  - `int` : return the solution and the derivatives up to that order,
    the parameter `return_derivatives=int(N)` is equivalent to `return_derivatives=range(N+1)`
  - `array_like` : specifies the list of derivative orders to be evaluated, 
    the 0-th derivative corresponds to the solution itself. 
    For example `return_derivatives=(0,1)` tells to return solution and it's first derivative
    and `return_derivatives=(1,)` tells to return only first derivative without solution itself
  - empty_array : return nothing, usefull when only return_delayed option is needed
* `return_delayed : None | "none" | "all" | "auto" | tuple(str)`
  parameter determines whether the delayed versions of solutions to be evaluated as well
  - `None | "none"` : do not return delayed variables
  - `"all"` : return all variables with all delays
  - `"auto"`: return only delayed variables that actually evaluated in the equation
  - `tuple(str)` : return specific delayed variables; each delayed 
    variable or its derivative is represented 
    by string `<variable>|<delay>`, `<variable>'|<delay>`, etc.
  note: the `return_derivatives` and `return_delayed` affect the output independently
* `stepsize : ("auto", dict) | "at" | array | float`
  - `("auto", dict)` : use adaptive stepsize control, dict specifies the parameters `tol`, `atol`, `rtol`, and `max_stepsize`:
    * `tol : float`
      sets `atol = rtol = tol` parameters for adaptive stepsize control
    * `atol : float`
      maximal estimated absolute error allowed at each step when using adaptive stepsize control
    * `rtol : float`
      maximal estimated relative error allowed at each step when using adaptive stepsize control
    * `max_stepsize : float`
      maximal allowed stepsize when using adaptive stepsize control
  - `"at"` : use the values specified by `at` keyword, valid only if `at` were specified by int or an array
  - `array` : use the values in array for integration points, 
    the parameters `stepsize=np.linspace(0,20,100), at="all"` are 
    equivalent to `stepsize='at', at=np.linspace(0,20,100)`
  - `float` : use fixed stepsize set by this value
* `max_steps : int`
  maximum number of steps the integrator is allowed to make, 
  if this number of steps is exceeded, `stop_integration` is 
  issued with a corresponding warning. Negative values or zero effectively disables this limitation.

* `discontinuities : str | tuple(str)`
  specifies the formulas when discontinuities may occur, str has the 
  form `"<lhs> = <rhs>"`, where `<lhs>` is formula containing variables and parameters,
  and `<rhs>` is does not contain any variables and parameters. Example: `discontinuities="x = 1"`.
  Note: discontinuities, introduced by `sign`, `abs`, `sign<a,b>`, `floor`, `Piecewise` functions 
  are handled automatically, for other cases use this option to include intersections with 
  discontinuity surfaces into integration mesh to avoid 
  order failure or excessive rejected steps in of numerical method.
* `discontinuity_detection` : bool
  default: true
* `events : oddesa.event | tuple(oddesa.event)`
* `method` : str
  numerical method to use, including rk4, rk45, and actually good methods




```python
oddesa.event(event : str, 
             condition : str or tuple(str),
             save : str or tuple(str),
             save_above : str or tuple(str),
             save_below : str or tuple(str),
             change : str or tuple(str),
             change_above : str or tuple(str),
             change_below : str or tuple(str),
             action : "step_on", "stop_integration" or "disable_event"
)
```

event could be string like "x = 0", or special like "step_rejected", "step"