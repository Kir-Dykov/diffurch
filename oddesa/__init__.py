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

class equation:

    def __repr__(self):
        return "equation(" + \
        "; ".join([f"{v}' = {rhs}" for v, rhs in self.equations.items()]) + \
        ")" + \
        f" where (" + \
        "; ".join([f"{p} = {val}" for p, val in self.parameters.items()]) + \
        ")"
    
    def __init__(self, 
                 equation : str, 
                 variational_equation : str = None, 
                 discontinuities : str = None, 
                 time_variable = "t", 
                 interval = None,
                 ic = None,
                 **kwargs):
        
        self.variable_names = list() # defined by the lhs of equations
        self.parameter_names = list() 
        
        self.equations = dict() # { var : it's derivative }
        
        self.delays = dict() # { delay: [affected variables] }
        self.delayed_variable_names = list()
        
        self.uuid = str(uuid.uuid4()).replace("-", "_")
        
        self.time_variable = time_variable
        
        self.highest_derivatives = list()
        
        self.interval = interval
        self.ic = ic
                
        # extacting lhs-variables and rhs
        for equality in equation.strip().split(";"):
            if not equality:
                continue
                
            lhs_rhs = equality.split("=")
            assert len(lhs_rhs) == 2, f"Wrong number of `=` signs in equation {repr(equality)}"
            lhs, rhs = lhs_rhs
            lhs_match = re.fullmatch(r"\s*([a-zA-Z_][a-zA-Z_0-9_]*'*)'\s*", lhs)
            assert lhs_match, f"The left-hand-side variable name in equation {repr(equality)} is incorrect."
            lhs_variable = lhs_match.group(1)
            self.variable_names.append(lhs_variable)
            self.highest_derivatives.append(lhs_variable+"'")
            self.equations[lhs_variable] = rhs.strip()
            # for higher order equations preform an order reduction, by using lower order derivatives as new variables
            while lhs_variable[-1] == "'":
                self.equations[lhs_variable[:-1]] = lhs_variable
                lhs_variable = lhs_variable[:-1]
                self.variable_names.append(lhs_variable) 
            
        # check if all equations are solved for highest order derivative (delayed derivatives are allowed in the right hand side: it is neutral type delay)
        for lhs, rhs in self.equations.items():
            for derivative in self.highest_derivatives:
                # x' = 2 x' (bad) # x' = 2 x'|1 (good)  # x'' = x''|1 (good) # x' = 2 x''|1 (bad)
                assert not re.search(f"{derivative}(?!\|)", rhs), f"Equation ({lhs}' = {rhs}) is not solved with respect to the highest derivatives."
        
        # detect parameter names as undefined tokens in gcc
        parameter_detection_code = "#include<math.h>\nint main() { \n\t" + \
            "\n\t".join([f"double {var};" for var in self.variable_names]) + "\n\t" + \
            "\n\t".join([f"double {var};" for var in self.highest_derivatives]) + "\n\t" + \
            "\n\t".join([f"{var} = {rhs};" for var, rhs in self.equations.items()]) + \
            "\n}"
        parameter_detection_code = parameter_detection_code.replace("'", "´").replace("|", "*")
        self.parameter_names = get_undefined_variables(parameter_detection_code)
        for parameter_name in self.parameter_names:
            assert parameter_name in kwargs, f"The parameter {repr(parameter_name)} is not specified."
            self.parameters[parameter_name] = kwargs[parameter_name]
        
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
        
        
        
        self.de_cpp_declaration = \
f"""struct MyDifferentialEquation : ??? {{
    {(chr(10)+4*" ").join([f"double {name};" for name, val in self.parameters.items()])}
    
    MyDifferentialEquation() {{
        f = [this](Vec<{len(self.variable_names)+len(self.delayed_variable_names)}> ∀){{
            const auto& [{", ".join(convert_identifiers(v) for v in self.variable_names + self.delayed_variable_names)}] = ∀;
            return Vec<{len(self.variable_names)}>{{
                {("," + chr(10)+16*" ").join([convert_identifiers(self.equations[var]) for var in self.variable_names])}
            }};
        }};
        
        delays = {{ 
            {{ {", ".join(f"Delay({delay})" for delay in self.delays.keys())} }},
        }}
    }}
}}
"""
    
    def print_cpp_parameters(self):
        cpp_code = f"""
        #include <iostream>
        using namespace std;
        int main() {{
            cout {" ".join([f' {chr(10)+16*" "}<< "{name}" << " = " << {val} << "; "' for name, val in self.parameters.items()])} << endl;
        }}
        """
        print(cpp_code)
        
        command = ['g++', '-Wall', '-Werror', '-x', 'c++', '-', '-o', 'print.o']
        result = subprocess.run(command, input=cpp_code, 
                            capture_output=True, text=True)
                          
        command = ['./print.o']
        result = subprocess.run(command, 
                            capture_output=True, text=True)
        
        print(result.stdout)
    
    def test001(self):
        cpp_code = f"""
        #include <iostream>
        #include "equation.hpp"
        using namespace std;
        int main() {{
            cout {" ".join([f' {chr(10)+16*" "}<< "{name}" << " = " << {val} << "; "' for name, val in self.parameters.items()])} << endl;
        
            auto eq = Lorenz(make_tuple(1,2,3));
            // eq.solution();
            cout << eq.solution() << endl;
            //printType(eq.solution());
        }}
        """
        print(cpp_code)
        
        command = ['pwd']
        result = subprocess.run(command, 
                                capture_output=True, text=True)
        print(result.stdout)
        
        command = ['g++', '-std=c++20', '-Wall', '-Werror', '-I', 'oddesa/cpp/', '-x', 'c++', '-', '-o', 'print.o']
        result = subprocess.run(command, input=cpp_code, 
                            capture_output=True, text=True)
        print(result.stderr)
                          
        command = ['./print.o']
        result = subprocess.run(command, 
                            capture_output=True, text=True)
        
        print(result.stdout)
        
    def test002(self):
        cpp_code = """
            #include <pybind11/pybind11.h>
            #include <pybind11/numpy.h>
            #include <vector>
            #include <tuple>

            namespace py = pybind11;

            // Function to separate even and odd numbers up to N
            std::tuple<py::array, py::array> separate_even_odd(int N) {
                std::vector<int> evens;
                std::vector<int> odds;

                for (int i = 1; i <= N; ++i) {
                    if (i % 2 == 0) {
                        evens.push_back(i);
                    } else {
                        odds.push_back(i);
                    }
                }

                // Convert std::vector to py::array
                py::array evens_array = py::array(evens.size(), evens.data());
                py::array odds_array = py::array(odds.size(), odds.data());

                return std::make_tuple(evens_array, odds_array);
            }

            // Binding code
            PYBIND11_MODULE(even_odd, m) {
                m.def("separate_even_odd", &separate_even_odd, "Separate numbers up to N into even and odd NumPy arrays");
            }
        """
        print(cpp_code)
        
        command = ['pwd']
        result = subprocess.run(command, 
                                capture_output=True, text=True)
        print(result.stdout)
        
        
        process = subprocess.run("g++ -O3 -Wall -Werror -shared \
        -std=c++20 -fPIC `python3 -m pybind11 --includes` \
        -o even_odd`python3-config --extension-suffix` \
        -x c++ -", capture_output=True, text=True, shell=True, input=cpp_code)
        
        print("STDERR:", process.stderr)
        print("STDOUT:", process.stdout)
        
        import even_odd
        evens, odds = even_odd.separate_even_odd(10)
        print("Evens:", evens)
        print("Odds:", odds)
        
        
        
        
    def compile(self):
        
        parameters_declaration = (chr(10)+16*" ").join(
            [f"double {name};" for name in self.parameters.keys()]
        )
        parameters_n = len(self.parameters.keys())
        parameters_list = ", ".join(self.parameters.keys())
        
        variables_n = len(self.variable_names)
        variables_list = ", ".join(self.variable_names)
        
        parameter_arguments_list = ", ".join([f"double {param}" for param in self.parameters.keys()])
        parameter_vector_arguments_list = ", ".join([f"const vector<double>& {param}" for param in self.parameters.keys()])
        
        cpp_code = f"""
            #include <pybind11/pybind11.h>
            #include <pybind11/stl.h>
            #include <pybind11/numpy.h>
            #include "equation.hpp"
            #include "utils/cartesian_product.hpp"
            #include <vector>
            #include <tuple>
            namespace py = pybind11;

            struct System : IVP<{variables_n}, System> {{
                {parameters_declaration}
                        
                System(const tuple<{", ".join(["double"]*parameters_n)}>& params) {{
                    tie({parameters_list}) = params;

                    lhs = [this](double t, const Vec<{variables_n}>& XXX){{
                        const auto& [{variables_list}] = XXX;
                        return Vec<{variables_n}>{{
                            {(","+chr(10)+28*" ").join(self.equations.values())}
                        }};
                    }};

                    initial_condition = 
                        [this](double t){{return Vec<{len(self.variable_names)}>{{{", ".join(f"{x}" for x in self.ic)}}};}};

                    integration_interval = make_pair({self.interval[0]}, {self.interval[1]});
                }}
                
                tuple<Event<tuple<{", ".join(["double"]*(variables_n+1))}>>> step_events = make_tuple(
                    Event<tuple<{", ".join(["double"]*(variables_n+1))}>> {{
                        .what_to_save = [this]() {{
                            const auto& [{variables_list}] = X;
                            return make_tuple(t, {variables_list});
                        }}
                    }}
                );
                
            }};
            
            auto solution_for_each({parameter_vector_arguments_list}) {{
                auto parameter_combinations = cartesian_product({parameters_list});
                
                vector<decltype(make_tuple(parameter_combinations[0], System(parameter_combinations[0]).solution()))> result;
                for (size_t i = 0; i < parameter_combinations.size(); i++) {{
                    auto& params = parameter_combinations[i];
                    auto eq = System(params);
                    result.push_back(make_tuple(params, eq.solution()));
                }}
                return result;
            }}
            
            auto solution({parameter_arguments_list}) {{
                auto params = make_tuple({parameters_list});
                auto eq = System(params);
                return make_tuple(params, eq.solution());
            }}

            // Binding code
            PYBIND11_MODULE(solution{self.uuid}, m) {{
                m.def("solution", &solution, 
                    "Calucate the solutions of the differential equaion."
                    {"".join([f', py::arg("{param}")' for param in self.parameters.keys()])});
            }}
        """
        cpp_code = convert_identifiers(cpp_code);
        
        print(cpp_code)
        
        command = ['pwd']
        result = subprocess.run(command, 
                                capture_output=True, text=True)
        print(result.stdout)
        
        
        process = subprocess.run(f"g++ -O3 -Wall -Werror -shared \
        -std=c++20 -fPIC  -Ioddesa/cpp/ `python3 -m pybind11 --includes` \
        -o solution{self.uuid}`python3-config --extension-suffix` \
        -x c++ -", capture_output=True, text=True, shell=True, input=cpp_code)
        
        print("STDERR:", process.stderr)
        print("STDOUT:", process.stdout)
        
        
    def solution(self, *args, **kwargs):
        
        solution_module = importlib.import_module(f'solution{self.uuid}')
        return solution_module.solution(*args, **kwargs)
        
        
        
        
        
        
        
        
        
        
#         syntax_error_detection_code = "#include<math.h>\nint main() { \n\t" + \
#             f"double {self.time_variable} = 0;\n\t" + \
#             "\n\t".join([f"auto {var} = " + "[=](double t){return double(0.);};" for var in self.variable_names]) + "\n\t" + \
#             "\n\t".join([f"double {par} = 0;" for par in self.parameter_names]) + "\n\t" + \
#             "\n\t".join([f"{rhs};" for var, rhs in self.equations.items()]) + \
#             "\n}"
        
#         errors = cpp_errors(syntax_error_detection_code)
#         if errors:
#             print(syntax_error_detection_code)
#             print(errors)
            
        """

        """