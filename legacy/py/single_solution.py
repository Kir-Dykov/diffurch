import json
import py.myutils as my
import importlib
importlib.reload(my)
import json
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.patheffects as path_effects


params = dict(
    b = 1, c = 1, d = -1, tau = 1, v = 1, 
    t_finish = 22, h = 0.01, comment = "real case, stable cycle"
)


params = dict(
    b = -0.2, c = 1, d = 1, tau = 1, v = 1,
    t_finish = 1000, h = 0.01, comment = "lorenz attractor"
)


params = dict(
    b = -0.02, c = 1, d = 1, tau = 1, v = 1,
    t_finish = 10, h = 0.01, comment = "restler attractor"
)


params = dict(
    b = -0.2, c = 1, d = 1, tau = 1, v = 1,
    t_finish = 50, h = 0.01, comment = "lorenz attractor"
)
#
script = "single_solution"
prefix = script
params_str = json.dumps(params)
my.run_cpp(script, params=params_str)
#
# just solution
filename = prefix + " " + params_str
x_, dx_, t_ = my.get_binary(f"output_bin/{filename}.bin")
plt.clf()
plt.plot(t_, x_)
plt.savefig(f"output_img/{filename}.jpg", format="jpg", dpi=1000)
os.system(f"kitten icat 'output_img/{filename}.jpg'")

# parametric
filename = prefix + " " + params_str + " (phase space projection)"
x_, dx_, t_ = my.get_binary(f"output_bin/{filename}.bin")
plt.clf()
plt.plot(x_, dx_)
plt.savefig(f"output_img/{filename}.jpg", format="jpg", dpi=1000)
os.system(f"kitten icat 'output_img/{filename}.jpg'")
