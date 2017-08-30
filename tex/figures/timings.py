import matplotlib.cm
import matplotlib.pyplot
import matplotlib.patches
import numpy
import pylatexenc.latexencode

# Global settings for matplotlib
matplotlib.pyplot.rc("text", usetex=True)
matplotlib.pyplot.rc("font", family="serif")
matplotlib.rcParams["hatch.linewidth"] = 0.5

cplex_x = [2, 4, 6, 8, 10]
cplex_y = [0.01145885, 0.02533634, 0.05837558, 0.12921851, 0.29617412]
cbc_x = [2, 4, 6, 8, 10]
cbc_y = [0.0225553491, 0.251695409, 1.76972775, 6.90875088, 41.8653300]
gplk_x = [2, 4, 6]
gplk_y = [0.00852113, 0.03818099, 1.25363812]

fig, ax = matplotlib.pyplot.subplots(figsize=(4,3))
ax.scatter(cplex_x, cplex_y, label="CPLEX", marker="^")
ax.scatter(cbc_x, cbc_y, label="Cbc", marker="o")
ax.scatter(gplk_x, gplk_y, label="GPLK", marker="*")
ax.set_yscale("log")
ax.set_xlabel("number of buffers, $N$")
ax.set_ylabel("solution time (seconds)")
ax.set_xticks(numpy.arange(0,11,2))
ax.set_ylim([0.001,100])
matplotlib.pyplot.legend()
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig("timings.pdf")
