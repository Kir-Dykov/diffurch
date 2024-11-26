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
from matplotlib.ticker import MaxNLocator



params = dict(
        size = 100,
        mode = "l1l2", 
        x_start = -5, x_finish = 5,
        y_start = -5, y_finish = 5,
        d = -1, tau = 1,
        v_start = 0.01, v_finish = 10, v_n = 50,
        t_finish = 25, h = 0.01,
)

recalculate = recompile = True

recalculate = recompile = False

params = dict(
        size = 50,
        mode = "reim", 
        x_start = -0.1, x_finish = 0.1,
        y_start = 0, y_finish = 4,
        d = -1, tau = 1,
        v_start = 0.01, v_finish = 30, v_n = 50,
        t_finish = 25, h = 0.01,
)

script = "fixed_points_count"
prefix = script
params_str = json.dumps(params)
#
filename = prefix + " " + params_str
filename_bin = f"output_bin/{filename}.bin"
#
if recalculate or not os.path.isfile(filename_bin):
    my.run_cpp(script, params=params_str, recompile=recompile)

# data
filename = prefix + " " + params_str
N = my.get_binary(f"output_bin/{filename}.bin").copy()
xs = np.linspace(params["x_start"], params["x_finish"], params["size"])
ys = np.linspace(params["y_start"], params["y_finish"], params["size"])
# figure
fig, ax = plt.subplots()
im = ax.pcolormesh(xs, ys, N, cmap="Grays")
plt.colorbar(im,ax=ax)
ax.set_aspect((params["x_finish"]-params["x_start"])/((params["y_finish"]-params["y_start"])))
# saving
fig.savefig(f"output_img/{filename}.jpg", format="jpg", dpi=1000, bbox_inches='tight', pad_inches=0.1)
os.system(f"kitten icat 'output_img/{filename}.jpg'")


re=-0.05
im=1.9
params.update(
        b = -2 * re,
        c = re*re + im*im
)
first_return_map(params, recalculate=recalculate, recompile=recompile)
