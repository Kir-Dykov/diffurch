import json
import os
import numpy as np
import shlex


output_path = "../output"
output_bin_path = output_path + "/bin"
cpp_path = "../cpp"
cpp_compiled_path = cpp_path + "/cpp_compiled"


# expected data format:
"""
number_of_arrays : int64
arrays : Array [number_of_arrays], where
    Array
        type_name : char, denotes the type of the array, the table is given in the code, type_size is decided by that table
        array_shape_size : int64, number of dimensions 
        array_shape : int64[array_shape_size], sizes of each dimensions
        array : type[...], flat array of type given by typename, and size given by array_shape

char: data_type (see table in the code)
"""
# TODO : string arrays support, like "Sc" to denote string of chars, and "Si" to denote arrays of ints of variable sizes. And each string is stored as [size, elem0, elem1, etc]
def get_binary(filename): 
    with open(f"{filename}", 'rb') as f:
        number_of_arrays = np.frombuffer(f.read(8), dtype='int64')[0]
        
        arrays = []
        
        for array_i in range(number_of_arrays):
        
            array_type = f.read(1).decode()
            type_name, type_size = {
                'f' : ('float32', 4),
                'd' : ('float64', 8), 
                'F' : ('float128', 16), 
                'c' : ('int8',    1),
                'i' : ('int32',   4),
                'u' : ('uint32',  4),
                'I' : ('int64',   8),
                'U' : ('uint64',  8),
            }[array_type]
            
        
            array_shape_size = np.frombuffer(f.read(8),                    dtype='int64')[0]
            array_shape      = np.frombuffer(f.read(8 * array_shape_size), dtype='int64')
            array            = np.frombuffer(f.read(type_size * np.prod(array_shape)), dtype=type_name).reshape(array_shape)
            arrays.append(array)
        
        return arrays

    
# def filename_from_params(params):
#     return json.dumps(params).replace('"',"").replace(": ","=")



# no .cpp extention in filename
def run_cpp(script_name, params, compiler_params = {}, has_output=True, recalculate=True, recompile=True, flags = "-O3"):
    
    params_str = shlex.quote(json.dumps(params))
    compiler_params_str =  ' '.join([f'-D{key}={shlex.quote(str(value))}' for key, value in compiler_params.items()])

    output_filename     = script_name + " " + params_str
    output_filename     = output_filename[:min(len(output_filename), 200)] # truncate by 200 characters
    output_filename_bin = f"{output_bin_path}/{output_filename}.bin"
    
    if recalculate or (has_output and not os.path.isfile(output_filename_bin)):    
        if not os.path.isdir(output_path):
            os.system(f"mkdir {output_path}")
        if not os.path.isdir(output_bin_path):
            os.system(f"mkdir {output_bin_path}")

        if recompile or not os.path.isfile(f"{cpp_compiled_path}/{script_name}.o"):
            if not os.path.isdir(cpp_compiled_path):
                os.system(f"mkdir {cpp_compiled_path}")

            cpp_compile_command = f"g++ -std=c++23 {flags} {compiler_params_str} -w {cpp_path}/{script_name}.cpp -o {cpp_compiled_path}/{script_name}.o"
            compile_status = os.system(cpp_compile_command)
            assert compile_status == 0, f"ERROR: compiler exited and returned {compile_status}"
        
        os.system(f'{cpp_compiled_path}/{script_name}.o {params_str} {shlex.quote(output_filename)}')
    
    if has_output:
        return get_binary(output_filename_bin)
    else:
        return None
                

