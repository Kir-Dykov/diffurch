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
#
#
#
def first_return_map(params, recalculate=True, recompile=True, scatter=True):
    script = "first_return_map"
    prefix = script
    params_str = json.dumps(params)
    #
    filename = prefix + " " + params_str
    filename_bin = f"output_bin/{filename}.bin"
    #
    if recalculate or not os.path.isfile(filename_bin):
        my.run_cpp(script, params=params_str, recompile=recompile)
    #
    # image
    v_, p_, t_return, zero_count = my.get_binary(filename_bin).copy()
    v_ = np.where(p_ == 0, np.nan, v_)
    # p_ = np.where(p_ == 0, np.nan, p_)
    # p_ = np.where(p_ == 0, np.nan, p_)
    diff_ = np.abs(np.diff(p_, prepend=p_[0]))
    discontinuities = (diff_ > 10*np.mean(diff_)) 
    p0 = np.where(discontinuities, np.nan, p_)
    t0 = np.where(discontinuities, np.nan, t_return)
    z0 = np.where(discontinuities, np.nan, zero_count)
    print(np.max(t_return))
    print(np.max(zero_count))
    t_return = np.where(t_return == 0, np.nan, t_return)
    fig, ax = plt.subplots()
    #
    ax.plot(v_, p_, c='b', lw=1, linestyle=':')
    ax.plot(v_, p0, c='b', lw=2)
    if scatter:
        ax.scatter(v_, p0, c='k', s=2)
    ax.plot(v_, v_*0, c='k')
    ax.plot(v_, v_,   c='k')
    ax.set_ylabel('$|p|$', color='b', rotation=0)
    ax.tick_params(axis='y', labelcolor='b')
    ax.dataLim.y1 = np.nanmax(p0)*1.05
    # ax.dataLim.y1 = 2 
    #
    ax_t = ax.twinx()
    ax_t.plot(v_, t_return, c='r', lw=1, linestyle=':')
    ax_t.plot(v_, t0, c='r', lw=2)
    ax_t.set_ylabel('t_0', color='r', rotation=0)
    ax_t.tick_params(axis='y', labelcolor='r')
    ax_t.dataLim.y0 = 0 
    #
    ax_z = ax.twinx()
    ax_z.plot(v_, z0, c='orange')
    ax_z.plot(v_, zero_count, c='orange', lw=1, linestyle=':')
    ax_z.set_ylabel('zero count', color='orange')
    ax_z.tick_params(axis='y', labelcolor='orange')
    ax_z.spines.right.set_position(("axes", 1.15))
    ax_z.dataLim.y0 = 0 
    ax_z.yaxis.set_major_locator(MaxNLocator(integer=True))
    fig.savefig(f"output_img/{filename}.jpg", format="jpg", dpi=1000, bbox_inches='tight', pad_inches=0.1)
    os.system(f"kitten icat 'output_img/{filename}.jpg'")
    plt.clf()




params = dict(
    b = 1, c = 1, d = -1, tau = 1, 
    v_start = 0.001, v_finish = 10, v_n = 100,
    t_finish = 22, h = 0.01, comment = "real case, stable cycle"
)

params = dict(
    b = -0.2, c = 1, d = 1, tau = 1,
    v_start = 0.001, v_finish = 3, v_n = 10000,
    t_finish = 1000, h = 0.01, comment = "lorenz attractor"
)


params = dict(
    b = -0.2, c = 1, d = 1, tau = 1,
    v_start = 0.5, v_finish = 0.7, v_n = 1000,
    t_finish = 100, h = 0.0001, comment = "lorenz attractor"
)


params = dict(
    b = -0.02, c = 1, d = 1, tau = 1,
    v_start = 0.001, v_finish = 6, v_n = 1000,
    t_finish = 100, h = 0.01, comment = "restler attractor"
)

params = dict(
    b = -0.02, c = 1, d = 1, tau = 1,
    v_start = 0.001, v_finish = 0.6, v_n = 1000,
    t_finish = 100, h = 0.001, comment = "restler attractor"
)



params = dict(
    b = -0.2, c = 1, d = 1, tau = 1,
    v_start = 0.001, v_finish = 3, v_n = 10000,
    t_finish = 100, h = 0.001, comment = "lorenz attractor output"
)

params = dict(
    b = 1, c = 1, d = -1, tau = 1, 
    v_start = 0.001, v_finish = 3, v_n = 100,
    t_finish = 22, h = 0.001, comment = "real case, stable cycle"
)

re = -0.25; im = 2
params = dict(
    b = -2*re, c = re*re +  im*im, d = -1, tau = 1, 
    v_start = 0.001, v_finish = 10, v_n = 100,
    t_finish = 22, h = 0.001, comment = "complex case"
)


recalculate = recompile = True

recalculate = recompile = False

params = dict(
    b = -0.2, c = 1, d = 1, tau = 1,
    v_start = 0.1, v_finish = 3, v_n = 100000,
    t_finish = 100, h = 0.01, comment = "lorenz attractor"
)

first_return_map(params, recalculate=recalculate, recompile=recompile)




#solutions
plt.clf()
fig, ax = plt.subplots()
for v in np.linspace(0.3,0.4,9):
    params["v"] = v
    params["h"] = 0.2
    params["t_finish"] = 10 
    #
    script = "single_solution"
    prefix = script
    params_str = json.dumps(params)
    my.run_cpp(script, params=params_str)
    #
    # just solution
    filename = prefix + " " + params_str
    x_, dx_, t_ = my.get_binary(f"output_bin/{filename}.bin")
    ax.plot(t_, x_)
    ax.scatter(t_, x_,s=4,c='k')
    ax.plot(t_, 0*t_, c='k') 
    \r\r\r
plt.savefig(f"output_img/{filename}.jpg", format="jpg", dpi=1000)
os.system(f"kitten icat 'output_img/{filename}.jpg'")

# parametric
filename = prefix + " " + params_str
x_, dx_, t_ = my.get_binary(f"output_bin/{filename}.bin")
filename += " (parametric)"
plt.clf()
plt.plot(x_, dx_)
plt.savefig(f"output_img/{filename}.jpg", format="jpg", dpi=1000)
os.system(f"kitten icat 'output_img/{filename}.jpg'")
script = "single_solution"
prefix = script
params_str = json.dumps(params)
my.run_cpp(script, params=params_str)
