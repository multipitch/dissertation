#!/usr/bin/env python

import argparse
import configparser
import csv
import os

import cplex
import numpy


parser = argparse.ArgumentParser(prog="model.py", 
                                 description="Run buffer vessel simulation.")
parser.add_argument("-P", "--path", default=".", type=str,
                    help="file path (default: <current directory>)")
parser.add_argument("-p", "--parameters", default="parameters.ini", type=str,
                    help="parameters filename (default: parameters.ini)")
parser.add_argument("-b", "--buffers", default="buffers.csv", type=str,
                    help="buffers filename (default: buffers.csv)")
parser.add_argument("-v", "--vessels", default="vessels.csv", type=str,
                    help="vessel filename (default: vessels.csv)")
parser.add_argument("-t", "--tolerance", default="1E-6", type=float,
                    help="tolerance for boolean values (default: 1e-06)")

# globals
ARGS = parser.parse_args()
PATH = os.path.abspath(os.path.expanduser(ARGS.path)) + "/"
PARAMETERS = PATH + ARGS.parameters
BUFFERS = PATH + ARGS.buffers
VESSELS = PATH + ARGS.vessels
TOLERANCE = ARGS.tolerance


class Parameters:
    
    def __init__(self, filename=PARAMETERS):
        file_parser = configparser.SafeConfigParser()
        file_parser.read(filename)
        parameters = file_parser.items("parameters")
        self.input_dict = dict([(p[0], float(p[1])) for p in parameters])
        self.cycle_time = self.input_dict["cycle_time"]
        self.prep_pre_duration = self.input_dict["prep_pre_duration"]
        self.prep_post_duration = self.input_dict["prep_post_duration"]
        self.transfer_duration = self.input_dict["transfer_duration"]
        self.hold_pre_duration = self.input_dict["hold_pre_duration"]
        self.hold_post_duration = self.input_dict["hold_post_duration"]
        self.hold_duration_min = self.input_dict["hold_duration_min"]
        self.hold_duration_max = self.input_dict["hold_duration_max"]
        self.minimum_fill_ratio = self.input_dict["minimum_fill_ratio"]
    
    def __repr__(self):
        lines = ""
        for k, v in self.__dict__.items():
            if k != "input_dict":
                lines += "{}: {}\n".format(k,v)
        return lines[:-1] # strip trailing newline


class Vessels:
    
    def __init__(self, filename=VESSELS):
        self.filename = filename
        self.input_dict = csv_columns_to_dict_of_lists(self.filename)
        self.names = self.input_dict["names"]
        self.volumes = self.input_dict["volumes"]
        self.costs = self.input_dict["costs"]
        self.max_volume = max(self.volumes)
        self.count = len(self.names)


class Buffers:
    
    def __init__(self, cycle_time, filename=BUFFERS):
        self.filename = filename
        self.input_dict = csv_columns_to_dict_of_lists(self.filename)
        self.names = self.input_dict["names"]
        self.volumes = self.input_dict["volumes"]
        self.absolute_use_start_times = self.input_dict["use_start_times"]
        self.use_durations = self.input_dict["use_durations"]
        self.relative_use_start_times = []
        self.use_start_times_offset_cycles = []
        for t in self.absolute_use_start_times:
            self.relative_use_start_times.append(t % cycle_time)
            self.use_start_times_offset_cycles.append(t // cycle_time)
        self.use_start_times = self.relative_use_start_times # shortened name
        self.count = len(self.names)

class Constraints:
    
    def __init__(self):
        self.rhs = []
        self.senses = []
        self.names = []
        self.coefficients = []
        self.counter = 0

    def add(self, rhs, senses, names, coefficients):
        self.rhs.append(rhs)
        self.senses.append(senses)
        self.names.append(names)
        for i in coefficients:
            self.coefficients.append((self.counter, *i))
        self.counter += 1

class Variables:
    
    def __init__(self):
        self.names = []
        self.obj = []
        self.lb = []
        self.ub = []
        self.types = []
        self.counts = []
        self.dimensions_list = []
        self.offsets = []
        self.short_names = []
        self.count = 0
    
    def add(self, name, dimensions, obj, lb, ub, types):
        if name in self.names:
            raise ValueError('name already exists')
        self.short_names.append(name)
        if len(self.offsets):
            self.offsets.append(self.offsets[-1] + self.counts[-1])
        else:
            self.offsets.append(0)
        count = prod(dimensions) 
        self.counts.append(count)
        self.dimensions_list.append(dimensions)
        self.names.extend(generate_names_from_dimensions(name, dimensions))
        for item in [obj, lb, ub, types]:
            if len(item) != count:
                raise ValueError('inconsistent data sizes')
        self.obj.extend(obj)
        self.lb.extend(lb)
        self.ub.extend(ub)
        self.types.extend(types)
        self.count += 1
    
    def pos(self, name, index):
        try:
            name_position = self.short_names.index(name)
        except ValueError:
            print ('name does not exist')
            raise
        offset = self.offsets[name_position]
        dimensions = self.dimensions_list[name_position]
        if type(index) == int and type(dimensions) == int:
            return offset + index
        len_index = len(index)
        if len_index != len(dimensions):
            raise ValueError('index is the wrong shape')
        for i, num in enumerate(index):
            if num not in range(dimensions[i]):
                raise IndexError('index is out of range')
        if type(index) == int:
            return offset + index
        multipliers = [1]
        for n in range(len_index - 1):
            multipliers.append(dimensions[len_index - 1 - n] * multipliers[n])
        multipliers = multipliers[::-1]
        rel_pos = sum([index[i] * multipliers[i] for i in range(len_index)])
        return offset + rel_pos


def csv_columns_to_dict_of_lists(filename):
    with open(filename) as f:
        reader = csv.reader(f, skipinitialspace=True, delimiter=",",
                            quoting=csv.QUOTE_NONNUMERIC)
        data_list = list(reader)
        lines = len(data_list)
        data_dict = {}
        for i, key in enumerate(data_list[0]):
            values = []
            for j in range(1, lines):
                values.append(data_list[j][i])
            data_dict[key] = values
        return data_dict


def prod(dimensions):
    try:
        product = 1
        for i in dimensions:
            product *= i
        return product
    except:
        return dimensions


def generate_names_from_dimensions(name, dimensions):
    # TODO: this can probably be simplified - see pos method in Variable
    names = []
    try:
        ndims = len(dimensions)
        periods = [1]
        count = 1
        for d in range(ndims):
            count *= dimensions[d]
            if d == ndims - 1:
                break
            periods.append(periods[-1] * dimensions[1 - d])
        periods = periods[::-1]
        for c in range(count):
            new_name = name + "_("
            for d in range(ndims):
                new_name += "{},".format((c // periods[d]) % dimensions[d])
            names.append(new_name[:-1] + ")")
    except:
        count = dimensions
        for n in range(dimensions):
            names.append("{}_{}".format(name, n))        
    return names


def build_variables(parameters, buffers, vessels):
    variables = Variables()
    M = vessels.count
    N = buffers.count
    P = N # number of slots
    
    # initialise booleans x_(n,p): buffer n made in slot p for all n in N for 
    # all p in P
    variables.add(name = "x",
                  dimensions = (N, P),
                  obj = [0] * N * P,
                  lb = [0] * N * P,
                  ub = [1] * N * P,
                  types = ["B"] * N * P)
    
    # initialise booleans y_(m,p): vessel size m in slot p for all m in M for
    # all p in P
    variables.add(name="y",
                  dimensions = (M, P),
                  obj = [vessels.costs[m] for m in range(M)] * P,
                  lb = [0] * M * P,
                  ub = [1] * M * P,
                  types = ["B"] * M * P)
    
    # initialise z_n: hold duration for buffer n for all n in N, incl the 
    # following bounds: 
    # parameters.hold_duration_min <= z_n <= parameters.hold_duration_max
    variables.add(name = "z",
                  dimensions = (N),
                  obj = [0] * N,
                  lb = [parameters.hold_duration_min] * N,
                  ub = [parameters.hold_duration_max] * N,
                  types = ["C"] * N)
    
    # initialise indicator booleans w_(n,k,p): buffers n and k both made in 
    # slot p for all n in N for all k in N for all p in P
    variables.add(name = "w",
                  dimensions = (N, N, P),
                  obj = [0] * N * N * P,
                  lb = [0] * N * N * P,
                  ub = [1] * N * N * P,
                  types = ["B"] * N * N * P)
    
    # initialise indicator booleans v_(n,k): buffers n and k share a slot
    # for all n in N for all k in N
    variables.add(name = "v",
                  dimensions = (N, N),
                  obj = [0] * N * N,
                  lb = [0] * N * N,
                  ub = [1] * N * N,
                  types = ["B"] * N * N)
    
    # initialise indicator booleans u_(n,k): buffer n is preparead after 
    # buffer k for all n in N for all k in N
    variables.add(name = "u",
                  dimensions = (N, N),
                  obj = [0] * N * N,
                  lb = [0] * N * N,
                  ub = [1] * N * N,
                  types = ["B"] * N * N)
    
    return variables


def build_constraints(parameters, buffers, vessels, variables):
    c =  Constraints()
    pos = variables.pos    
    M = vessels.count
    N = buffers.count
    P = N  # number of slots
    
    # each buffer can only be made in one slot
    for n in range(N):
        coeffs = [(pos("x", (n, p)), 1) for p in range(P)]
        c.add(rhs=1, 
              senses="E",
              names="one_slot_per_buffer_{}".format(n),
              coefficients=coeffs)
    
    # each slot contains a max of one vessel
    for p in range(P):
        coeffs = [(pos("y", (m, p)), 1) for m in range(M)]
        c.add(rhs=1,
              senses="L",
              names="max_one_vessel_per_slot_{}".format(p),
              coefficients=coeffs)
    
    # vessel must be large enough
    for n in range(N):
        for p in range(P):
            coeffs = [(pos("x", (n, p)), buffers.volumes[n])]
            for m in range(M):
                coeffs.append((pos("y", (m, p)), 0 - vessels.volumes[m]))
            c.add(rhs=0,
                  senses="L",
                  names="slot_{}_vessel_big_enough_for_buffer_{}".format(p,n),
                  coefficients=coeffs)
    
    # vessel must be small enough
    for n in range(N):
        for p in range(P):
            coeffs = [(pos("x",(n,p)), vessels.max_volume)]
            for m in range(M):
                coeffs.append((pos("y", (m, p)),
                        vessels.volumes[m] * parameters.minimum_fill_ratio))
            c.add(rhs=buffers.volumes[n] + vessels.max_volume,
                  senses="L",
                  names="slot_{}_vessel_small_enough_for_buffer_{}".format(p,n),
                  coefficients=coeffs)
    
    # total duration in buffer hold must not exceed cycle time
    for n in range(N):
        rhs = (parameters.cycle_time - parameters.hold_pre_duration
               - parameters.hold_post_duration - buffers.use_durations[n])
        c.add(rhs=rhs,
              senses="L",
              names="buffer_{}_hold_vessel_duration_<=_cycle_time".format(n),
              coefficients=[(pos("z", n), 1)])
    
    # track if buffers n and k both made in slot p
    for n in range(N):
        for k in range(N):
            if k != n:
                for p in range(P):
                    coeffs = [(pos("x", (n, p)), 1), 
                              (pos("x", (k, p)), 1),
                              (pos("w", (n, k, p)), -1)]
                    names="buffers_{}_and_{}_share_slot_{}".format(n, k, p)
                    c.add(rhs=1,
                          senses="L",
                          names=names,
                          coefficients=coeffs)
    
    # track if buffers n and k both made in same slot
    for n in range(N):
        for k in range(n + 1, N):
            coeffs = [(pos("v", (n, k)), -1)]
            for p in range(P):
                coeffs.append((pos("w", (n, k, p)), 1))
            c.add(rhs=0,
                  senses="L",
                  names="buffers_{}_and_{}_share_a_slot".format(n, k),
                  coefficients=coeffs)
    
    # prep scheduling - only one prep can take place in a slot at a given time
    prep_duration = (parameters.prep_pre_duration 
                     + parameters.transfer_duration
                     + parameters.prep_post_duration)
    for n in range(N):
        for k in range(n + 1, N):
            for i in (0, 1):
                rhs = ((2 * i - 1) * buffers.use_start_times[n]
                       + (1 - 2 * i) * buffers.use_start_times[k]
                       + 2 * (1 + i) * parameters.cycle_time
                       - prep_duration)
                names = ("scheduling_buffer_{}_after_buffer_{}"
                         .format((k, n)[i], (n, k)[i]))
                coeffs = [(pos("z", n), 2 * i - 1),
                          (pos("z", k), 1 - 2 * i),
                          (pos("u", (n, k)), 
                           2 * (2 * i - 1) * parameters.cycle_time),
                          (pos("v", (n, k)), 2 * parameters.cycle_time)]
                c.add(rhs=rhs, senses="L", names=names, coefficients=coeffs)
    
    return c


def add_variables(prob, variables):
    prob.variables.add(variables.obj, variables.lb, variables.ub, 
                       variables.types, variables.names)


def add_constraints(prob, linear_constraints):
    prob.linear_constraints.add(rhs=linear_constraints.rhs,
                                senses=linear_constraints.senses,
                                names=linear_constraints.names)
    prob.linear_constraints.set_coefficients(linear_constraints.coefficients)


# check buffer data TODO: this is currently unused!!!
def check_buffer_data():    
    if len(buffers.use_start_times) != N or len(buffers.use_durations) != N:
        raise ValueError ("Buffer data incomplete")
    for n in buffers.use_start_times:
        if n >= parameters.cycle_time:
            raise ValueError("use_start_times must be modulo cycle time")
    for n in buffers.volumes:
        if n > max(vessels.volumes):
            raise ValueError("buffer volume(s) greater than max vessel volume")


def solve_model(prob):
    try:
        prob.solve()
    except cplex.exceptions.CplexSolverError as e:
        print("Exception raised during solve: " + e)
    else:
        solution = prob.solution
        print("\nSolution status = ", solution.get_status(), ":", end=' ')
        print(solution.status[solution.get_status()])


def print_solution_variables():
    vnames = prob.variables.get_names()
    for i, j in enumerate(vnames):
        print("{}\t:  {}".format(j, prob.solution.get_values(i)))


# TODO: get rid of this ASAP - use numpy or pandas to present results nicely
def print_solution(prob,variables):
    pos = variables.pos
    values = prob.solution.get_values()
    selected_slots = []
    selected_volumes = []
    for m in range(vessels.count):
        for p in range(buffers.count):
            position = pos("y", (m,p))
            if abs(values[position] - 1.0) <= TOLERANCE:
                selected_slots.append(p)
                selected_volumes.append(vessels.volumes[m])
    lines = ["|slots    |"]
    for s in selected_slots: lines[0] += " {}\t\t|".format(s)
    lines.append("|volumes  |")
    for s in selected_volumes: lines[1] += " {}\t|".format(s)
    lines.append("\n|buffers  ")
    for n in range(buffers.count):
        line = "| {}\t  |".format(n)
        for p in selected_slots:
            position = pos("x", (n,p)) 
            if abs(values[position] - 1.0) <= TOLERANCE:
                line += " X\t\t|"
            else:
                line += "  \t\t|"
        lines.append(line)
    for l in lines: print(l)


# TODO: reimplement as a Results class
def unflatten_results(parameters, buffers, vessels, constraints, solutions):
    M = vessels.count
    N = buffers.count
    P = N  # number of slots
    results = {}
    for index, variable in enumerate(variables.short_names):
        start = variables.offsets[index]
        end = start + variables.counts[index]
        results[variable] = numpy.reshape(numpy.asarray(solutions[start:end]),
                                          variables.dimensions_list[index])
        if variables.types[index] =="B":
            results[variable] = numpy.rint(results[variable]).astype(int)
    return results


def simple_plot(parameters, buffers, vessels, constraints, results):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    # TODO: add prep vessels, label all vessels by volume and type, 
    # add legend for buffers, with name and prep volume
    
    M = vessels.count
    N = buffers.count
    P = N  # number of slots
    colors=list(cm.Set2(numpy.linspace(0,1,N)))
    vessels, slots = numpy.nonzero(results["y"])
    
    fig, ax = plt.subplots()
    for n in range(N):
        hold_start_time = (buffers.use_start_times[n] 
                           - results["z"][n]
                           - parameters.transfer_duration
                           - parameters.hold_pre_duration)
        hold_duration = (parameters.hold_pre_duration
                         + parameters.transfer_duration
                         + results["z"][n]
                         + buffers.use_durations[n]
                         + parameters.hold_post_duration)
        if hold_start_time < 0:
            xranges = [(0, hold_start_time + hold_duration),
                       (hold_start_time + parameters.cycle_time,
                        parameters.cycle_time)]
        else:
            xranges = [(hold_start_time, hold_duration)] 
        
        bar_height = 0.6
        ax.broken_barh(xranges, (N - (0.5 + n + 0.5 * bar_height), bar_height),
                       facecolors=colors[n], edgecolors="black", zorder=3)
    #ax.grid(True, zorder=0)        
    ax.grid(axis="x", linestyle="solid", linewidth=1, zorder=0)
    ax.grid(axis="y", linestyle="dashed", linewidth=1, zorder=0)
    ax.set_ylim(0, N)
    ax.set_xlim(0, parameters.cycle_time)
    ax.set_xlabel('time (h)')
    ax.set_ylabel('Buffer Hold Vessels')
    ax.set_yticks([n + 0.5 for n in range(N)])
    ax.set_xticks([6 * (t + 1) for t in range(int(parameters.cycle_time / 6))])
    hold_names = buffers.names
    ax.set_yticklabels([name for name in buffers.names][::-1])
    plt.show()

if __name__ == "__main__":
    parameters = Parameters()
    buffers = Buffers(parameters.cycle_time)
    vessels = Vessels()

    variables = build_variables(parameters, buffers, vessels)
    constraints = build_constraints(parameters, buffers, vessels, variables)    
    
    prob = cplex.Cplex()
    prob.objective.set_sense(prob.objective.sense.minimize)
    
    add_variables(prob, variables)
    add_constraints(prob, constraints)
    
    solve_model(prob)
    
    solutions = prob.solution.get_values()
    results = unflatten_results(parameters, buffers, vessels, constraints,
                                  solutions)
    simple_plot(parameters, buffers, vessels, constraints, results)
    
