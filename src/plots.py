import matplotlib.cm
import matplotlib.pyplot
import numpy
import scipy.stats


def single_cycle_plot(parameters, buffers, vessels, filename=None):    
    # TODO: add legend, naming convention changes, include hatching / sub-bars
    # TODO: function for generating list of (start, duration) pairs based on
    # an input (start, duration) and cycle time
    # TODO: make more use of numpy / vectorized functions to remove some loops
    # TODO: sorting/ranking below is messy - should be able to simplify
    
    M = vessels.count
    N = buffers.count
    P = N  # number of slots    
    ct = parameters.cycle_time    
    colors=list(matplotlib.cm.tab20(numpy.linspace(0,1,N)))    
    used_slots = set(buffers.prep_slots)
    nslots = len(used_slots)
    sorted_slots = sorted(buffers.p_to_m, key = buffers.p_to_m.get)
    slot_ranks = {}
    for j, i in enumerate(sorted_slots):
        slot_ranks[i] = j
    sorted_prep_names = [buffers.p_to_m[i] for i in sorted_slots]
    prep_index = [slot_ranks[i] for i in buffers.prep_slots] 
    bar_height = 0.6
    fig, ax = matplotlib.pyplot.subplots(figsize=(10, (N + nslots + 1) / 4))
    
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
                           edgecolors="black", linewidth=1, zorder=3)
    
    # Use Bars
    for n in range(N):
        xranges = cyclic_xranges(buffers.use_start_times[n],
                                 buffers.use_durations[n], ct)
        ystart =  N +  nslots - (prep_index[n] - 0.5 + 0.5 * bar_height)
        yrange = (N - (0.5 + n + 0.5 * bar_height), bar_height)
        ax.broken_barh(xranges, yrange, facecolors=colors[n], hatch="\\\\\\",
                       edgecolors="black", linewidth=1, zorder=3)
    
    # Procedure Outlines
    for n in range(N):
        ax.broken_barh(hold_xranges[n], hold_yranges[n], facecolors='none',
                       edgecolors="black", linewidth=1, zorder=4)
        ax.broken_barh(prep_xranges[n], prep_yranges[n], facecolors='none',
                       edgecolors="black", linewidth=1, zorder=4)
    
    # Axes
    ax.grid(axis="x", linestyle="solid", linewidth=1, zorder=0)
    ax.grid(axis="y", linestyle="dashed", linewidth=1, zorder=0)
    ax.set_xlabel('time (h)')
    ax.set_ylabel('Vessels')
    ax.set_yticks([n + 0.5 for n in range(N)]
                  + [i + 0.5 for i in range (N + 1, N + nslots + 2)])
    ax.set_xticks([6 * (t + 1) for t in range(int(parameters.cycle_time / 6))])
    prep_lbl = ["{} Prep ".format(vessels.names[i]) for i in sorted_prep_names]
    hold_lbl = ["{} Hold".format(n) for n in buffers.names]
    ax.set_yticklabels(hold_lbl[::-1] + prep_lbl[::-1])
    ax.set_xlim(0, parameters.cycle_time)
    ax.set_ylim(0, N + nslots + 1)

    # Title
    matplotlib.pyplot.title("Equipment Time Utilisation for One Cycle")
    
    # Write to file or plot to screen
    matplotlib.pyplot.tight_layout()
    if filename:
        matplotlib.pyplot.savefig(filename)
    else:
        matplotlib.pyplot.show()


# If cyclic operation continues past end of cycle, split into two operations
def cyclic_xranges(start_time, duration, ct):
    if 0 <= start_time < ct:
        if start_time + duration > ct:
            return [(0, duration - (ct - start_time)), 
                    (start_time, ct - start_time)]
        return [(start_time, duration)]
    raise ValueError("start time out of range")
