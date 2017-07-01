#!/usr/bin/env python

import argparse
import configparser
import csv
import os

import numpy
import pulp


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
parser.add_argument("-s", "--solver", default="", type=str,
                    help="solver to be used (default: <pulp default solver>")


# globals
ARGS = parser.parse_args()
PATH = os.path.abspath(os.path.expanduser(ARGS.path)) + "/"
PARAMETERS = PATH + ARGS.parameters
BUFFERS = PATH + ARGS.buffers
VESSELS = PATH + ARGS.vessels
SOLVER = ARGS.solver


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


class Variable:
    
    def __init__(self, name, dimensions, low_bound, up_bound, cat):
        self.name = name
        self.dimensions = dimensions
        self.low_bound = low_bound
        self.up_bound = up_bound
        self.cat = cat
        self.objects = unflatten_variable(
                pulp.LpVariable.dicts(self.name,
                                      tuple([range(i) for i in dimensions]),
                                      self.low_bound, self.up_bound, self.cat),
                self.dimensions)
        self.values = numpy.empty(self.dimensions)
        self.o = self.objects # shortened name
    
    def evaluate(self):
        vectorized_evaluation = numpy.vectorize(lambda i: pulp.value(i))
        self.values = vectorized_evaluation(self.objects)
        return self.values


class Variables:
    
    def __init__(self, varbunch=None):
        if varbunch:
            self.add(varbunch)
    
    def add(self, varbunch):
        for i in varbunch:
            self.__dict__[i.name] = i


class Results:
    
    def __init__(self, variables=None):
        if variables:
            self.add(variables)
    
    def add(self, variables):
        for name, var in variables.__dict__.items():
            self.__dict__[name] = var.evaluate()


# Generate a dict of lists from a csv file
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


# Return the product of the elements of a list, tuple or set
def prod(dimensions):
    try:
        product = 1
        for i in dimensions:
            product *= i
        return product
    except:
        return dimensions


# Generator for accessing components of a variable
def variable_iterator(variable, dimensions):
    if len(dimensions) > 1:
        for i in range(dimensions[0]):
            yield from variable_iterator(variable[i], dimensions[1:])
    else:
        for i in range(dimensions[0]):
            yield variable[i]


# Generator for accessing components of a variable, returning the component and
# its position
def variable_iterator_loc(variable, dimensions):
    index = [0]* len(dimensions)
    depth = 0
    yield from _variable_iterator_loc(variable, dimensions, index, depth)


# Recursive function called by variable_iter_loc
def _variable_iterator_loc(variable, dimensions, index, depth):
    depth += 1
    if len(dimensions) == 1:
        for i in range(dimensions[0]):
            index[depth - 1] = i
            yield variable[i], tuple(index)
    else:
        for i in range(dimensions[0]):
            index[depth - 1] = i
            yield from _variable_iterator_loc(variable[i], dimensions[1:],
                                              index, depth)


# create a k-dimensional numpy array of variable objects, given the variable
# as a (dict of dict of ... of dicts) of depth k.
def unflatten_variable(variable, dimensions):
    unflattened = numpy.empty(dimensions, object)
    variter = variable_iterator_loc(variable, dimensions)
    for i in range(prod(dimensions)):
        varobject, position = next(variter)
        unflattened[position] = varobject
    return unflattened


def define_problem(parameters, buffers, vessels):
    M = vessels.count
    N = buffers.count
    P = N  # number of slots
    
    ms = list(range(M))
    ns = list(range(N))
    ps = list(range(P))
    
    # Problem
    problem = pulp.LpProblem("Buffer Preparation Vessel Selection",
                             pulp.LpMinimize)
    
    # Variables
    u = Variable("u", (N, N), 0, 1, "Integer")
    v = Variable("v", (N, N), 0, 1, "Integer")
    w = Variable("w", (N, N, P), 0, 1, "Integer")
    x = Variable("x", (N, P), 0, 1, "Integer")
    y = Variable("y", (M, P), 0, 1, "Integer")
    z = Variable("z", (N,), parameters.hold_duration_min,
                 parameters.hold_duration_max, "Continuous")
    
    # Keep track of variables in a class
    variables = Variables([u, v, w, x, y, z])
    
    # Objective Function
    problem += sum([vessels.costs[m] * y.o[m][p] for m in ms for p in ps])
    
    # Constraint 1: A buffer may only be prepared in one slot
    for n in ns:
        problem += sum([x.o[n][p] for p in ps]) == 1
    
    # Constraint 2: Each slot can contain at most one vessel
    for p in ps:
        problem += sum([y.o[m][p] for m in ms]) <= 1
    
    # Constraint 3: Vessels must be big enough for buffers
    for n in ns:
        for p in ps:
            problem += (buffers.volumes[n] * x.o[n][p]
                        <= sum([vessels.volumes[m] * y.o[m][p] for m in ms]))
    
    # Constraint 4: Vessels must be small emough for buffers
    mfr = parameters.minimum_fill_ratio
    for n in ns:
        for p in ps:
            problem += (vessels.max_volume * x.o[n][p]
                        + mfr * sum([vessels.volumes[m]*y.o[m][p] for m in ms])
                        <= buffers.volumes[n]
                        + vessels.max_volume)
    
    # Constraint 5: Total hold proceudre durations must be less than cycle time
    for n in ns:
        rhs = (parameters.cycle_time
               - parameters.hold_pre_duration
               - parameters.hold_post_duration
               - buffers.use_durations[n])
        problem += z.o[n] <= rhs
    
    # Constraint 6: For each pair of buffers, indicate if made in same slot
    for n in ns:
        for k in range(n + 1, N):
            for p in ps:
                problem += x.o[n][p] + x.o[k][p] - w.o[k][n][p] <= 1
    
    # Constraint 7: For each pair of buffers, indicate which prepared first
    for n in ns:
        for k in range(n + 1, N):
            problem += sum([w.o[n][k][p] for p in ps]) - v.o[n][k] <= 0
    
    # Constraint 8: Each prep vessel can only do one thing at a time
    for i in (0, 1):
        for n in ns:
            for k in range(n + 1, N):
                problem += ((2 * i - 1) * z.o[n]
                            - (2 * i - 1) * z.o[k]
                            + (2 * i - 1) * 2*parameters.cycle_time * u.o[n][k]
                            + 2 * parameters.cycle_time * v.o[n][k]
                            - (2 * i - 1) * buffers.use_start_times[n]
                            + (2 * i - 1) * buffers.use_start_times[k]
                            - (i + 1) * 2 * parameters.cycle_time
                            + parameters.prep_pre_duration
                            + parameters.transfer_duration
                            + parameters.prep_post_duration
                            <= 0)
    
    return problem, variables


if __name__ == "__main__":
    parameters = Parameters()
    buffers = Buffers(parameters.cycle_time)
    vessels = Vessels()
    problem, variables = define_problem(parameters, buffers, vessels)
    
    problem.writeLP(PATH + "model.lp")
    
    if SOLVER == "CPLEX":
        problem.solve(pulp.solvers.CPLEX())
    else:
        problem.solve()
    
    print("Status:", pulp.LpStatus[problem.status])
    
    results = Results(variables)
    
    from plots import single_cycle_plot
    single_cycle_plot(parameters, buffers, vessels, results,
                      filename=(PATH + "plot.pdf"))
