import json
import os
import numpy as np

# expected data format:
# byte number | meaning
# 0 | array type: single char
#   | f - float32, F - float64 (double),  
#   | i - int32,   I - int64 (long)
#   | u - uint32,  U - uint64
#   | c - int8 (char)
# +4 | number of dimensions: int32
# +4*number of dimensions | shape: int32*
# + all data flattened
def get_binary(filename): 
    with open(f"{filename}", 'rb') as f:
        array_type = f.read(1).decode()
        type_name, type_size = {
            'f' : ('float32', 4),
            'd' : ('float64', 8),
            'i' : ('int32', 4),
            # 'I' : ('int64', 8),
            # 'u' : ('uint32', 4),
            # 'U' : ('uint64', 8),
            # 'c' : ('int8', 1),
        }[array_type]
        
        array_shape_size = np.frombuffer(f.read(4), dtype='int32')[0]
        # print(array_shape_size)
        array_shape = np.frombuffer(f.read(4 * array_shape_size), dtype='int32')
        # print(array_shape)
        array = np.frombuffer(f.read(type_size * np.prod(array_shape)), dtype=type_name).reshape(array_shape)
        
        return array

    
def filename_from_params(params):
    return json.dumps(params).replace('"',"").replace(": ","=")

# no .cpp extention in filename
def run_cpp(filename, params = "{}", output_prefix = None, recompile=True):
    if output_prefix is None:
        output_prefix = filename

    if recompile or not os.path.isfile(f"cpp_compiled/{filename}.o"):
        cpp_compile_command = f"g++ -std=c++20 -Wfatal-errors -w cpp/{filename}.cpp -o cpp_compiled/{filename}.o"
        exit_code = os.system(cpp_compile_command)
        if exit_code == 0: 
            os.system(f'./cpp_compiled/{filename}.o \'{params}\' \'{output_prefix}\'')
        else:
            print(f"ERROR: compiler exited and returned {exit_code}")
    else:
        os.system(f'./cpp_compiled/{filename}.o \'{params}\' \'{output_prefix}\'')

