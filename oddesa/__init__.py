import subprocess
import re
import regex
import importlib
import uuid


acute_accent = "´" # for derivative variable name
dental_click = "ǀ" # Latin Letter Dental Click ǀ # will be used to denote delayed argument
ogonek = "˛" + "ˌ"

unicode_token_replacements = [
    ("'", "´"), 
    ("|", "ǂ"), 
    # (".", "‸"),
]

def convert_identifiers(a: str):
    for FROM, TO in unicode_token_replacements:
        a = a.replace(FROM, TO)
    return a

def get_undefined_variables(source_code):
    command = ['g++', '-Wall', '-Werror', 
               '-fsyntax-only', '-x', 'c++', '-']
    result = subprocess.run(command, input=source_code, 
                            capture_output=True, text=True)
    pattern = r"error: ‘([^‘]+)’ was not declared in this scope"
    undefined_vars = re.findall(pattern, result.stderr)
    # print(result.stderr)
    return undefined_vars

def cpp_errors(source_code):
    command = ['g++', '-Wall', '-Werror', 
               '-Wno-unused-value', '-Wno-unused-variable', 
               '-fsyntax-only', '-x', 'c++', '-']
    result = subprocess.run(command, input=source_code, 
                            capture_output=True, text=True)
    return result.stderr

class event:
    def __init__(self, where, filter=None, save=None, change=None, action=None):
        self.where = where
        self.filter = filter
        self.save = save
        self.change = change
        self.action = action
        
    def __repr__(self):
        return f"event{self.__dict__}"

class equation:
    def __repr__(self):
        return f"equation{self.__dict__}"
        
    def __init__(self, 
                 equation : str, 
                 at = "all",
                 events = None,
                 # variational_equation : str = None, 
                 # discontinuities : str = None, 
                 # independent_variable = "t", 
                 compile = True,
                 **kwargs):
        
        # self.independent_variable = independent_variable 
        self.equations = dict() # { var : it's derivative }
        self.initial_conditions = dict() # { var : formula for initial condition }
        
        self.variable_names = list()
        self.parameter_names = list() 
        self.delays = dict() # { delay: [affected variables] }
        self.delayed_variable_names = list()
        self.highest_derivatives = list()
        
        
        # Universally unique identifier, used for compiled package name
        # to avoid confusion, and imposible-to-reimport issue
        self.uuid = ""
                
        # extacting lhs-variables and rhs
        for equality in equation.strip().split(";"):
            if not equality:
                continue
                
            lhs_rhs = equality.split("=")
            assert len(lhs_rhs) == 2, f"Wrong number of `=` signs in equation {repr(equality.strip())}"
            lhs, rhs = lhs_rhs
            rhs = rhs.strip()
            lhs_match = re.fullmatch(r"\s*([a-zA-Z_][a-zA-Z_0-9_]*)('+)\s*", lhs) # x''' -> group(1) == "x", group(2) == "'''"
            assert lhs_match, f"The left-hand-side variable name in equation {repr(equality)} is incorrect."
            lhs_variable, primes = lhs_match.group(1), lhs_match.group(2)
            assert lhs_variable not in self.variable_names, f"Repeated definition of variable {lhs_variable} in equality {repr(equality.strip())}."
            
            self.highest_derivatives.append(lhs_variable + primes)
            
            for derivative_order in range(1, len(primes)):
                var = lhs_variable+"'"*(derivative_order-1)               
                self.variable_names.append(var)
                self.equations[var] = var+"'"
            
            var = lhs_variable + primes[:-1]
            self.variable_names.append(var)
            self.equations[var] = rhs
            

        # check if all equations are solved for highest order derivative (delayed derivatives are allowed in the right hand side: it is neutral type delay)
        for lhs, rhs in self.equations.items():
            for derivative in self.highest_derivatives:
                # x' = 2 x' (bad) # x' = 2 x'|1 (good)  # x'' = x''|1 (good) # x' = 2 x''|1 (bad)
                assert not re.search(f"{derivative}(?!\|)", rhs), f"System is not solved with respect to the highest derivatives, in equality ({lhs}' = {rhs})."
        
        # initial conditions handling
        assert 'initial_condition' not in kwargs and 'ic' not in kwargs, f"Duplicate initial condition specification using both keyword arguments 'initial_condition' and 'ic'."
        if 'initial_condition' in kwargs or 'ic' in kwargs:
            ic_string = kwargs['ic'] if 'ic' in kwargs else kwargs['initial_condition']
        else:            
            # implicit definition of initial condition as constant
            ic_string = "; ".join(f"{variable} = {variable}₀" for variable in self.variable_names) 
            
        for initial_condition in ic_string.strip().split(";"):
            if not initial_condition:
                continue
            lhs_rhs = initial_condition.split("=")
            assert len(lhs_rhs) == 2, f"Wrong number of `=` signs in initial condition {repr(initial_condition.strip())}"
            lhs, rhs = lhs_rhs
            var = lhs.strip()
            rhs = rhs.strip()
            assert var in self.variable_names, f"Ininital condition for an unknown variable {repr(var)} in {initial_condition.strip()}."
            self.initial_conditions[var] = rhs            
            
        for var in self.variable_names:
            assert var in self.initial_conditions, f"Initial condition for variable {repr(var)} is not specified."
        
        # detect parameter names as undefined tokens in gcc
        parameter_detection_code = "#include<math.h>\nint main() { \n\t" + \
            "\n\t".join([f"double {var};" for var in self.variable_names]) + "\n\t" + \
            "\n\t".join([f"double {var};" for var in self.highest_derivatives]) + "\n\t" + \
            "\n\t".join([f"{var} = {rhs};" for var, rhs in self.equations.items()]) + \
            "\n\t".join([f"{var} = {rhs};" for var, rhs in self.initial_conditions.items()]) + \
            "\n}"
        parameter_detection_code = parameter_detection_code.replace("'", "´").replace("|", "*")
        self.parameter_names = get_undefined_variables(parameter_detection_code)
                
        # for parameter_name in self.parameter_names:
        #     assert parameter_name in kwargs, f"The parameter {repr(parameter_name)} is not specified."
        #     self.parameters[parameter_name] = kwargs[parameter_name]
        
        # now search for delayed arguments
        pattern = r"([a-zA-Z_d][a-zA-Z_0-9_]*'*)\|([a-zA-Z_0-9_\.]+)"  # pattern = r'|\|\(((?:[^()]++|\((?1)\))*)\)'
        for lhs, rhs in self.equations.items():
            for var, delay in regex.findall(pattern, rhs):
                if delay not in self.delays:
                    self.delays[delay] = [var]
                elif var not in self.delays[delay]:
                    self.delays[delay].append(var)
                
                delayed_variable = f"{var}|{delay}"
                if delayed_variable not in self.delayed_variable_names:
                    self.delayed_variable_names.append(delayed_variable)
                    
        if events is None:
            self.events = []
        elif isinstance(events, list):
            self.events = events
        else:
            self.events = [events]
            
        match at:
            case "all":
                self.events.append(event("step", save=", ".join(["t"] + self.variable_names)))
            case "last":
                self.events.append(event("stop_integration", save=", ".join(["t"] + self.variable_names)))
            case None:
                pass
            case "none":
                pass
            case _:
                raise ValueError("Invalid value for the keyword 'at'.")
        
        print(self.events)
        if compile:
            self.compile(**kwargs)   
            
    def __del__(self):
        # delete shared object (.so) file created by self.compile
        if self.uuid:
            subprocess.run(f"rm solution{self.uuid}`python3-config --extension-suffix`", 
                           capture_output=True, text=True, shell=True)
        
    def compile(self, debug = False):
        
        if self.uuid:
            subprocess.run(f"rm solution{self.uuid}`python3-config --extension-suffix`", 
                           capture_output=True, text=True, shell=True)
        
        self.uuid = str(uuid.uuid4()).replace("-", "_")
        
        
        parameter_declaration = (chr(10)+16*" ").join(
            [f"double {name};" for name in self.parameter_names]
        )
        parameter_n = len(self.parameter_names)
        parameter_list = ", ".join(self.parameter_names)
        
        variable_n = len(self.variable_names)
        variable_list = ", ".join(self.variable_names)
        
        parameter_arguments_list = ", ".join([f"double {param}" for param in ["initial_time", "final_time"] + self.parameter_names])
        parameter_vector_arguments_list = ", ".join([f"const vector<double>& {param}" for param in self.parameter_names])
        
        
        
        cpp_code = \
f"""
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "equation.hpp"
#include "vector_to_numpy.hpp"
#include "utils/cartesian_product.hpp"
#include <vector>
#include <tuple>
namespace py = pybind11;

struct System : IVP<{variable_n}, System> {{
    {parameter_declaration}

    System(const tuple<{", ".join(["double"]*parameter_n)}>& params) {{
        tie({parameter_list}) = params;

        lhs = [this](double t, const Vec<{variable_n}>& XXX){{
            const auto& [{variable_list}] = XXX;
            return Vec<{variable_n}>{{
                {(","+chr(10)+28*" ").join(self.equations.values())}
            }};
        }};

        initial_condition = 
            [this](double t){{return Vec<{len(self.variable_names)}>{{{", ".join(f"{self.initial_conditions[var]}" for var in self.variable_names)}}};}};

    }}

    tuple<Event<tuple<{", ".join(["double"]*(variable_n+1))}>>> step_events = make_tuple(
        Event<tuple<{", ".join(["double"]*(variable_n+1))}>> {{
            .what_to_save = [this]() {{
                const auto& [{variable_list}] = X;
                return make_tuple(t, {variable_list});
            }}
        }}
    );

    // delays = {{ 
    //    {{ {", ".join(f"Delay({delay})" for delay in self.delays.keys())} }},
    // }}

}};

auto solution_for_each({parameter_vector_arguments_list}) {{
    auto parameter_combinations = cartesian_product({parameter_list});

    vector<decltype(make_tuple(parameter_combinations[0], System(parameter_combinations[0]).solution(0., 1.)))> result;
    for (size_t i = 0; i < parameter_combinations.size(); i++) {{
        auto& params = parameter_combinations[i];
        auto eq = System(params);
        result.push_back(make_tuple(params, eq.solution(0, 1)));
    }}
    return result;
}}

auto solution({parameter_arguments_list}) {{
    auto params = make_tuple({parameter_list});
    auto eq = System(params);
    return convert_to_numpy(make_tuple(params, eq.solution(initial_time, final_time)));
}}

// Binding code
PYBIND11_MODULE(solution{self.uuid}, m) {{
    m.def("solution", &solution, 
        "Calucate the solutions of the differential equaion."
        {"".join([f', py::arg("{param}")' for param in ["initial_time", "final_time"] + self.parameter_names])});
}}
"""
        cpp_code = convert_identifiers(cpp_code);
        
        process = subprocess.run(f"g++ -O3 -Wall -Werror -shared \
        -std=c++20 -fPIC  -Ioddesa/cpp/ `python3 -m pybind11 --includes` \
        -o solution{self.uuid}`python3-config --extension-suffix` \
        -x c++ -", capture_output=True, text=True, shell=True, input=cpp_code)
        
        if process.stderr:
            
            print("Compiler returned some errors.")
            print("Code:")
            
            for i, line in enumerate(cpp_code.split("\n")):
                print(i+1, line)
            
            print("STDERR of compiler:\n", process.stderr)
        
        
        
    def solution(self, interval, **kwargs):
        
        if not self.uuid:
            self.compile()
            
        solution_module = importlib.import_module(f'solution{self.uuid}')
        
        kwargs_ = dict()
        for key, value in kwargs.items():
            if convert_identifiers(key) in self.parameter_names:
                kwargs_[convert_identifiers(key)] = value
            if convert_identifiers(key) in self.variable_names:
                if isinstance(value, tuple):
                    for i, val in enumerate(value):
                        kwargs_[convert_identifiers(key + "'"*i + "₀")] = val
                else:
                    kwargs_[convert_identifiers(key + "₀")] = value
        # print(kwargs_)
        
        for param in self.parameter_names:
            assert param in kwargs_, f"Parameter {param} were not specified."
            
        
        return solution_module.solution(interval[0], interval[1], **kwargs_)
        
        
        