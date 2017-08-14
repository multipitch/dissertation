import csv

import matplotlib.cm
import matplotlib.pyplot
import numpy
import pylatexenc.latexencode
from scipy.misc import comb

# Global settings for matplotlib
matplotlib.pyplot.rc("text", usetex=True)
matplotlib.pyplot.rc("font", family="serif")
matplotlib.rcParams["hatch.linewidth"] = 0.5


def single_cycle_plot(parameters, buffers, vessels, filename=None):    

    M = vessels.count
    N = buffers.count
    P = N  # number of slots    
    ct = parameters.cycle_time
       
    colors = list(matplotlib.cm.tab20(numpy.linspace(0,1,N)))
    
    used_slots = set(buffers.prep_slots)
    nslots = len(used_slots)
    sorted_slots = sorted(buffers.p_to_m, key = buffers.p_to_m.get)
    slot_ranks = {}
    for j, i in enumerate(sorted_slots):
        slot_ranks[i] = j
    sorted_prep_names = [buffers.p_to_m[i] for i in sorted_slots]
    prep_index = [slot_ranks[i] for i in buffers.prep_slots] 
    bar_height = 0.6
    fig, ax = matplotlib.pyplot.subplots(figsize=(10, (N + nslots + 2) / 2.5))
    
    # Buffer Hold Vessel Bars
    hold_xranges = []
    hold_yranges = []
    for n in range(N):
        hold_xranges.append((cyclic_xranges(buffers.hold_start_times[n], 
                                            buffers.hold_total_durations[n],
                                            ct)))
        hold_yranges.append((N - (0.5 + n + 0.5 * bar_height), bar_height))               
    for n in range(N):
        ax.broken_barh(hold_xranges[n], hold_yranges[n], facecolors=colors[n],
                       zorder=3)
    
    # Buffer Prep Vessel Bars
    prep_xranges = []
    prep_yranges = []
    for n in range(N):
        prep_xranges.append(cyclic_xranges(buffers.prep_start_times[n],
                                           buffers.prep_total_durations[n],
                                           ct))
        ystart =  N + nslots - (prep_index[n] - 0.5 + 0.5 * bar_height)
        prep_yranges.append((ystart, bar_height))
    for n in range(N):
        ax.broken_barh(prep_xranges[n], prep_yranges[n], facecolors=colors[n],
                       zorder=3)
    
    # Tx Bars
    for n in range(N):
        t_tx = parameters.transfer_duration
        xranges = cyclic_xranges(buffers.transfer_start_times[n], t_tx, ct)
        ystart =  N + nslots - (prep_index[n] - 0.5 + 0.5 * bar_height)
        yranges = [(ystart, bar_height), 
                   (N - (0.5 + n + 0.5 * bar_height), bar_height)]
        for yrange in yranges:
            ax.broken_barh(xranges, yrange, facecolors=colors[n], hatch="///",
                           edgecolors="black", linewidth=0.5, zorder=3)
    
    # Use Bars
    for n in range(N):
        xranges = cyclic_xranges(buffers.use_start_times[n],
                                 buffers.use_durations[n], ct)
        ystart =  N +  nslots - (prep_index[n] - 0.5 + 0.5 * bar_height)
        yrange = (N - (0.5 + n + 0.5 * bar_height), bar_height)
        ax.broken_barh(xranges, yrange, facecolors=colors[n], hatch="\\\\\\",
                       edgecolors="black", linewidth=0.5, zorder=3)
    
    # Procedure Outlines
    for n in range(N):
        ax.broken_barh(hold_xranges[n], hold_yranges[n], facecolors="none",
                       edgecolors="black", linewidth=0.75, zorder=4)
        ax.broken_barh(prep_xranges[n], prep_yranges[n], facecolors="none",
                       edgecolors="black", linewidth=0.75, zorder=4)
    
    # Axes and Labels
    ax.grid(axis="x", linestyle="solid", linewidth=1, zorder=0)
    ax.grid(axis="y", linestyle="dashed", linewidth=1, zorder=0)
    ax.set_xlabel("time (h)")
    ax.set_ylabel("Vessels")
    ax.set_yticks([n + 0.5 for n in range(N)]
                  + [i + 0.5 for i in range (N + 1, N + nslots + 2)])
    ax.set_xticks([6 * (t + 1) for t in range(int(parameters.cycle_time / 6))])
    prep_labels = []
    for i in sorted_prep_names:
        prep_label = "{} Prep ".format(vessels.names[i])
        prep_labels.append(pylatexenc.latexencode.utf8tolatex(prep_label))
    hold_labels = []
    for n in buffers.names:
        hold_label = "{} Hold".format(n)
        hold_labels.append(pylatexenc.latexencode.utf8tolatex(hold_label))
    ax.set_yticklabels(hold_labels[::-1] + prep_labels[::-1])
    ax.set_xlim(0, parameters.cycle_time)
    ax.set_ylim(0, N + nslots + 1)
    matplotlib.pyplot.title("Steady-State Equipment Time Utilisation")
    
    # Write to file or plot to screen
    matplotlib.pyplot.tight_layout()
    if filename:
        matplotlib.pyplot.savefig(filename)
        matplotlib.pyplot.close("all")
    else:
        matplotlib.pyplot.show()


def explanatory_plot():
    
    fig, ax = matplotlib.pyplot.subplots(figsize=(10, 3))
    
    ax.broken_barh([(6, 10), (28,13)], (150, 20), facecolors="white",
                   edgecolors="black", linewidth=0.5, zorder=3)
    ax.broken_barh([(16, 12)], (90, 80), facecolors="white",
                   edgecolors="black", linewidth=0.5, zorder=3, 
                   linestyle="dotted")
    ax.broken_barh([(16, 12)], (150, 20), facecolors="white",
                   edgecolors="black", linewidth=0.5, zorder=3, hatch="///")
    ax.broken_barh([(6, 35)], (150, 20), facecolors="none",
                   edgecolors="black", linewidth=0.75, zorder=4)
    ax.broken_barh([(5, 11), (28, 36), (84, 11)], (90, 20),
                   facecolors="white", edgecolors="black", linewidth=0.5,
                   zorder=3)
    ax.broken_barh([(16, 12)], (90, 20), facecolors="white",
                   edgecolors="black", linewidth=0.5, zorder=3, hatch="///")
    ax.broken_barh([(64, 20)], (30, 80), facecolors="white",
                   edgecolors="black", linewidth=0.5, zorder=3, 
                   linestyle="dotted")
    ax.broken_barh([(64, 20)], (90, 20), facecolors="white",
                   edgecolors="black", linewidth=0.5, zorder=3, hatch="\\\\\\")
    ax.broken_barh([(5, 90)], (90, 20), facecolors="none",
                   edgecolors="black", linewidth=0.75, zorder=4)
    ax.broken_barh([(-10, 120)], (30, 20), facecolors="white",
                   edgecolors="none", linewidth=0.75, zorder=3) 
    ax.broken_barh([(64, 3), (69, 2), (73, 5), (80, 4)], (30, 20),
                   facecolors="white", edgecolors="black", linewidth=0.5,
                   zorder=3, hatch="\\\\\\")  
    ax.broken_barh([(-10, 120)], (30, 20), facecolors="none",
                   edgecolors="black", linewidth=0.75, zorder=4,
                   linestyle="dashed") 
                       
    ax.set_ylim(0, 200)
    ax.set_xlim(0, 100)
    ax.grid(axis="x", linestyle="solid", linewidth=1, zorder=0)
    ax.grid(axis="y", linestyle="dashed", linewidth=1, zorder=0)
    ax.set_xlabel("time (h)")
    ax.set_ylabel("Process Equipment")
    ax.set_yticks([40, 100, 160])
    ax.set_xticks([0, 64, 100])
    ax.set_xticklabels(["0", "$t_{USE,n}$", "$\lambda$"])
    ax.set_yticklabels(["(Process Users)", "Hold Vessel", 
                        "Preparation Vessel"])
                        
    ax.text(11, 160, "$\Delta t_{\mathit{PREP,PRE}}$", fontsize=8,
            horizontalalignment="center", verticalalignment="center")
    ax.text(34.5, 160, "$\Delta t_{\mathit{PREP,POST}}$", fontsize=8,
            horizontalalignment="center", verticalalignment="center")
    ax.text(10.5, 100, "$\Delta t_{\mathit{HOLD,PRE}}$", fontsize=8,
            horizontalalignment="center", verticalalignment="center")
    ax.text(89.5, 100, "$\Delta t_{\mathit{HOLD,POST}}$", fontsize=8,
            horizontalalignment="center", verticalalignment="center")
    ax.text(22, 130, "$\Delta t_{\mathit{TRANSFER}}$", fontsize=8,
            horizontalalignment="center", verticalalignment="center")
    ax.text(74, 70, "$\Delta t_{\mathit{USE},n}$", fontsize=8,
            horizontalalignment="center", verticalalignment="center")
    ax.text(48, 100, "$\\textbf{\\textit{z}}_{n}$", fontsize=8,
            horizontalalignment="center", verticalalignment="center")
            
    matplotlib.pyplot.title("Steady-State Equipment Time Utilisation")
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig("explanatory.pdf")
    matplotlib.pyplot.close("all")

def sched_plot_single():
    
    reftime = 60
    duration = 20
    offset = -15
    overhang = 10
    xmin = 0
    xmax = 100
    ymin = 0
    ymax = 100
    thickness = 20
    
    fig, ax = matplotlib.pyplot.subplots(figsize=(8, 2))
    
    ymid = 0.5 * (ymin + ymax)
    bottom = ymid - 0.5 * thickness    
    leftstart = xmin - overhang
    leftduration = reftime - duration - leftstart
    leftmid = 0.5 * (xmin + reftime - duration)
    rightstart = reftime + duration
    rightduration = xmax + overhang - rightstart
    rightmid = 0.5 * (reftime + duration + xmax)
    prepmid = reftime + offset + 0.5 * duration
    tprepn = "$\\textrm{\\boldmath$t$}_{\mathit{PREP},n}$"
    tprepk = "$\\textrm{\\boldmath$t$}_{\mathit{PREP},k}$"
    tupper = tprepn[:-1] + " - \Delta t_{\mathit{PREP}}$"
    tlower = tprepn[:-1] + " + \Delta t_{\mathit{PREP}}$"
    ax.broken_barh([(leftstart, leftduration), (rightstart, rightduration)],
                   (bottom, thickness), facecolors="0.75", linestyle="None",
                   zorder=3)
    ax.broken_barh([(reftime + offset, duration)], (bottom, thickness),
                   facecolors="white", edgecolors="black", linewidth=0.75,
                   zorder=3)
    ax.plot([reftime], [bottom + thickness], '.k', zorder=4)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin, xmax)
    ax.grid(axis="x", linestyle="solid", linewidth=1, zorder=0)
    ax.grid(axis="y", linestyle="dashed", linewidth=1, zorder=0)
    ax.set_xlabel("time (h)")
    ax.set_yticks([ymid])
    ax.set_xticks([xmin, reftime - duration, reftime, reftime + duration,
                   xmax])
    ax.set_xticklabels(["0", tupper, tprepn, tlower, "$T$"])               
    ax.set_yticklabels(["vessel $p$"])
    ax.text(leftmid, ymid, tprepk + " feasible region", fontsize=8,
            horizontalalignment="center", verticalalignment="center", zorder=4)
    ax.text(prepmid, ymid, "buffer $n$ prep. procedure", fontsize=8,
            horizontalalignment="center", verticalalignment="center", zorder=4)
    ax.text(rightmid, ymid, tprepk + " feasible region", fontsize=8,
            horizontalalignment="center", verticalalignment="center", zorder=4)
    
            
    matplotlib.pyplot.title("Steady-State Equipment Time Utilisation")
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig("sched1.pdf")
    matplotlib.pyplot.close("all")

def sched_plot_all():
    
    ct = 100 # cycle time
    duration = 20
    offset = -15
    overhang = 10
    xmin = 0
    xmax = ct
    ymin = 0
    thickness = 15
    spacing = 50
    reftimes = [0,
                0.5 * duration,
                duration,
                0.5 * ct,
                ct - duration,
                ct - 0.5 * duration,
                ct]
    reftimes.reverse()
    fig, ax = matplotlib.pyplot.subplots(figsize=(8, 5))
    lines = len(reftimes)
    ymax = lines * spacing
    ymids = [(i + 0.5) * spacing for i in range(lines)]
    bottoms= [i - 0.5 * thickness for i in ymids]
    leftstart = xmin - overhang
    leftdurations = [i - duration - leftstart for i in reftimes]
    rightstarts = [i + duration for i in reftimes]
    rightdurations = [xmax + overhang - i for i in rightstarts]
    tprepn = "$\\textrm{\\boldmath$t$}_{\mathit{PREP},n}$"
    tupper = tprepn[:-1] + " - \Delta t_{\mathit{PREP}}$"
    tlower = tprepn[:-1] + " + \Delta t_{\mathit{PREP}}$"
    
    prep_xranges = []
    free_xranges = []
    for i in range(lines):
        prep_xranges.append(cyclic_xranges((reftimes[i] + offset) % ct,
                                           duration, ct))
        free_xranges.append(cyclic_xranges((reftimes[i] + duration) % ct,
                                           ct - 2 * duration, ct))
                                           
    
    for i in range(lines):
        ax.broken_barh(prep_xranges[i],
                       (bottoms[i], thickness), facecolors="white",
                       edgecolors="black", linewidth=0.75, zorder=3)

    for i in range(lines):
        ax.broken_barh(free_xranges[i],
                       (bottoms[i], thickness), facecolors="0.75",
                       linestyle="None", zorder=3)
        ax.plot([reftimes[i]], [bottoms[i] + thickness], '.k', zorder=4)
    
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin, xmax)
    ax.grid(axis="x", linestyle="solid", linewidth=1, zorder=0)
    ax.grid(axis="y", linestyle="dashed", linewidth=1, zorder=0)
    ax.set_xlabel("time (h)")
    ax.set_yticks(ymids)
    ax.set_xticks([xmin, xmax])
    ax.set_xticklabels(["0", "$T$"])
    ylabels = [tprepn[:-1] + " = 0$",
               "$0 < " + tprepn[1:-1] + " < " + tupper[1:],
               tprepn[:-1] + " = " + tupper[1:],
               tupper[:-1] + " < " + tprepn[1:-1] + " < " + tlower[1:],
               tprepn[:-1] + " = " + tlower[1:],
               tlower[:-1] + " < " + tprepn[1:-1] + " < T$",
               tprepn[:-1] + "= T$"]    
    ylabels.reverse()
    ax.set_yticklabels(ylabels)            
    matplotlib.pyplot.title("Steady-State Equipment Time Utilisation")
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig("sched2.pdf")
    matplotlib.pyplot.close("all")
    

def complexity_plot(parameter):
    
    fig, ax = matplotlib.pyplot.subplots(figsize=(5.5, 3))
    N = numpy.arange(1, 41, 1)
    
    # baselines
    ax.plot(N, N, "--k", linewidth=1, label="$N$")
    ax.scatter(N[1:], comb(N[1:], 2), label="${{N}\choose{2}}$", marker=".",
               color="k", s=2)
    ax.plot(N, N**2, "-.k", linewidth=1, label="$N^2$")
    ax.plot(N, N**3, "-k", dashes=[4, 2, 1, 2, 1, 2], linewidth=1, 
            label="$N^3$")
    
    if parameter == "dimensions":
        # basic model
        ax.scatter(N, 2 * N**2, label="basic model", marker="^")
        # complete model
        complete_dims = (N[1:] + 2) * comb(N[1:], 2) + 2 * (N[1:]**2 + N[1:])
        ax.scatter(N[1:], complete_dims, label="complete model", marker="o")
    elif parameter == "equations":
        # basic model
        ax.scatter(N, 2 * N**2 + 3 * N + 1, label="basic model", marker="^")
        # complete model
        complete_eqns = (3 * N[1:]**2  + 7 * comb(N[1:], 2) + 6 * N[1:] + 1)
        ax.scatter(N[1:], complete_eqns, label="complete model", marker="o")
    else:
        raise ValueError("{} is not a valid parameter".format(parameter))
    
    ax.set_yscale("log")
    ax.set_xlabel("number of buffers, $N$")
    ax.set_ylabel("number of {}".format(parameter))
    matplotlib.pyplot.legend(loc=0, ncol=3)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig("{}.pdf".format(parameter))
    matplotlib.pyplot.close("all")


def timing_plot(inputdata):
    if type(inputdata) is str:
        sizes, durations = read_durations(inputdata)
    else:
        (sizes, durations) = inputdata
    fig, ax = matplotlib.pyplot.subplots(figsize=(3, 3))
    ax.boxplot(durations.transpose(), labels=sizes)
    ax.set_yscale("log")
    ax.set_xlabel("number of buffers, $N$")
    ax.set_ylabel("solution time (seconds)")
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig("timing.pdf")
    matplotlib.pyplot.close("all")    
    

# If cyclic operation continues past end of cycle, split into two operations
def cyclic_xranges(start_time, duration, ct):
    if 0 <= start_time < ct:
        if start_time + duration > ct:
            return [(0, duration - (ct - start_time)), 
                    (start_time, ct - start_time)]
        return [(start_time, duration)]
    raise ValueError("start time out of range")


def read_durations(filename="durations.csv"):
    with open(filename) as f:
        reader = csv.reader(f)
        sizes = next(reader)
    durations = numpy.genfromtxt("durations.csv", delimiter=",", 
                                 skip_header=True)
    return sizes, durations


if __name__ == "__main__":
    
    
    explanatory_plot()
    sched_plot_single()
    sched_plot_all()
    complexity_plot("dimensions")
    complexity_plot("equations")
    timing_plot("durations.csv")
    
