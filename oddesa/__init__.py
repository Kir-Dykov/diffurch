import subprocess
import re
import regex

acute_accent = "´" # for derivative variable name
dental_click = "ǀ" # Latin Letter Dental Click ǀ # will be used to denote delayed argument
ogonek = "˛" + "ˌ"

unicode_token_replacements = [("'", "´"), ("|", "ǂ"), (".", "‸")]

def make_identifier(a: str):
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
                 **kwargs):
        self.variable_names = list() # defined by the lhs of equations
        self.parameter_names = list() 
        self.parameters = dict() # { name : value(s) }
        self.equations = dict() # { var : it's derivative }
        self.delays = dict() # { delay: [affected variables] }
        self.delayed_variable_names = list()
        
        self.time_variable = time_variable
        
        self.highest_derivatives = list()
        
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
        
        
        
        de_cpp_declaration = \
f"""struct MyDifferentialEquation : ??? {{
    {(chr(10)+4*" ").join([f"double {name};" for name, val in self.parameters.items()])}
    
    MyDifferentialEquation() {{
        f = [this](Vec<{len(self.variable_names)+len(self.delayed_variable_names)}> ∀){{
            const auto& [{", ".join(make_identifier(v) for v in self.variable_names + self.delayed_variable_names)}] = ∀;
            return Vec<{len(self.variable_names)}>{{
                {("," + chr(10)+16*" ").join([make_identifier(self.equations[var]) for var in self.variable_names])}
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
            
        