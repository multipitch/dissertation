import matplotlib.pyplot as plt
import numpy
from matplotlib import cm


def single_cycle_plot(parameters, buffers, vessels, constraints, results,
                      to_file=False, filename=None):
    
    # TODO: add legend, naming convention changes, include hatching / sub-bars
    # TODO: function for generating list of (start, duration) pairs based on
    # an input (start, duration) and cycle time
    
    M = vessels.count
    N = buffers.count
    P = N  # number of slots
    
    x = results["x"]
    y = results["y"]
    z = results["z"]
    
    ct = parameters.cycle_time
    
    colors=list(cm.Set2(numpy.linspace(0,1,N)))
    
    _, buffer_slots = numpy.nonzero(x)
    selected_vessels, selected_slots = numpy.nonzero(y)
    used_slots = len(selected_slots)
    slots_to_vessels = {}
    for i in range(used_slots):
        slots_to_vessels[selected_slots[i]] = selected_vessels[i]
    buffer_vessels = [slots_to_vessels[j] for j in buffer_slots]
    s = sorted(buffer_vessels)    
    prep_index = [s.index(k) for k in buffer_vessels]
    
    bar_height = 0.6
    fig, ax = plt.subplots(figsize=(14, 4))
    prep_duration = (parameters.prep_pre_duration
                     + parameters.transfer_duration
                     + parameters.prep_post_duration)
    
    def cyclic_xranges(start_time, duration, ct):
        if start_time < 0:
            if start_time + duration <= 0:
                xranges = [(start_time + ct, duration)]
            else:
                xranges = [(0, start_time + duration),
                           (start_time + ct, ct)]
        else:
            xranges = [(start_time, duration)]
        return xranges
    
    # Buffer Hold Vessel Bars
    hold_xranges = []
    hold_yranges = []
    for n in range(N):
        hold_start_times = (buffers.use_start_times[n] - z[n]
                            - parameters.transfer_duration
                            - parameters.hold_pre_duration)
        hold_durations = (parameters.hold_pre_duration
                          + parameters.transfer_duration + z[n]
                          + buffers.use_durations[n]
                          + parameters.hold_post_duration)
        hold_xranges.append((cyclic_xranges(hold_start_times, 
                                            hold_durations, ct)))
        hold_yranges.append((N - (0.5 + n + 0.5 * bar_height),
                            bar_height))               
    for n in range(N):
        ax.broken_barh(hold_xranges[n], hold_yranges[n], facecolors=colors[n], 
                       zorder=3)
    
    # Buffer Prep Vessel Bars
    prep_xranges = []
    prep_yranges = []
    for n in range(N):
        prep_start_time = (buffers.use_start_times[n]
                           - results["z"][n]
                           - parameters.transfer_duration
                           - parameters.prep_pre_duration)
        prep_xranges.append(cyclic_xranges(prep_start_time, prep_duration, ct))
        ystart =  N + used_slots + 1 - (0.5 + prep_index[n] + 0.5 * bar_height)
        prep_yranges.append((ystart, bar_height))
    for n in range(N):
        ax.broken_barh(prep_xranges[n], prep_yranges[n], facecolors=colors[n], 
                       zorder=3)
    
    # Tx Bars
    for n in range(N):
        dt = parameters.transfer_duration
        tx_start_time = (buffers.use_start_times[n] - z[n] - dt)
        xranges = cyclic_xranges(tx_start_time, dt, ct)
        ystart =  N + used_slots + 1 - (0.5 + prep_index[n] + 0.5 * bar_height)
        yranges = [(ystart, bar_height), 
                   (N - (0.5 + n + 0.5 * bar_height), bar_height)]
        for yrange in yranges:
            ax.broken_barh(xranges, yrange, facecolors=colors[n], hatch="////",
                           edgecolors="black", linewidth=1, zorder=3)
    
    # Use Bars
    for n in range(N):
        xranges = cyclic_xranges(buffers.use_start_times[n],
                                 buffers.use_durations[n], ct)
        ystart =  N +  used_slots - (0.5 + prep_index[n] + 0.5 * bar_height)
        yrange = (N - (0.5 + n + 0.5 * bar_height), bar_height)
        ax.broken_barh(xranges, yrange, facecolors=colors[n], hatch="\\\\\\",
                       edgecolors="black", linewidth=1, zorder=3)
    
    # All Bars outlines
    for n in range(N):
        ax.broken_barh(hold_xranges[n], hold_yranges[n], facecolors='none',
                       edgecolors="black", linewidth=1, zorder=4)
        ax.broken_barh(prep_xranges[n], prep_yranges[n], facecolors='none',
                       edgecolors="black", linewidth=1, zorder=4)


    
    ax.grid(axis="x", linestyle="solid", linewidth=1, zorder=0)
    ax.grid(axis="y", linestyle="dashed", linewidth=1, zorder=0)
    #ax.set_ylim(0, N + used_slots)
    ax.set_xlim(0, parameters.cycle_time)
    ax.set_xlabel('time (h)')
    ax.set_ylabel('Vessels')
    ax.set_yticks([n + 0.5 for n in range(N)]
                  + [i + 0.5 for i in range (N + 1, N + used_slots + 2)])
    ax.set_xticks([6 * (t + 1) for t in range(int(parameters.cycle_time / 6))])
    prep_lbls = ["{} Prep ".format(vessels.names[i]) for i in selected_vessels]
    hold_labels = ["{} Hold".format(n) for n in buffers.names]
    ax.set_yticklabels(hold_labels[::-1] + prep_lbls[::-1])
    plt.title("Equipment Time Utilisation for One Cycle")
    
    if filename:
        plt.savefig("plot.svg")
        #plt.savefig(filename)
    else:
        plt.show()
