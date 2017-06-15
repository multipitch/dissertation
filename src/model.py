import csv
import math
from configparser import SafeConfigParser

import cplex

import numpy as np


# global parameters
TOLERANCE = 1E-6


class Parameters:
    
    def __init__(self, filename="parameters.ini"):
        parser = SafeConfigParser()
        parser.read(filename)
        parameters = parser.items("parameters")
        self.__dict__ = dict([(p[0], float(p[1])) for p in parameters])
    
    def __repr__(self):
        lines = ""
        for k, v in self.__dict__.items():
            lines += "{}: {}\n".format(k,v)
        return lines[:-1] # strip trailing newline


class DataClass:
    
    def __init__(self, filename):
        self.filename = filename
        self.__dict__ = csv_columns_to_dict_of_lists(self.filename)
    
    def __repr__(self):
        lines = ""
        data_columns = []
        for k, v in self.__dict__.items():
            lines += "{:<20}".format(k)
            data_columns.append(v)
        rows = len(data_columns[0])
        columns = len(self.__dict__)
        if rows > 0:
            for r in range(rows):
                line = "\n"
                for c in range(columns):
                    line += "{:<20}".format(data_columns[c][r])
                lines += line
        return lines


class Vessels(DataClass):
    
    def __init__(self, filename="vessels.csv"):
        self.filename = filename
        self.__dict__ = csv_columns_to_dict_of_lists(self.filename)


class Buffers(DataClass):
    
    def __init__(self, filename="buffers.csv"):
        self.filename = filename
        self.__dict__ = csv_columns_to_dict_of_lists(self.filename)


class LinearConstraints:
    
    def __init__(self):
        self.rhs = []
        self.senses = []
        self.names = []
        self.coefficients = []


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
        
    def add(self, name, dimensions, obj, lb, ub, types):
        if name in self.names:
            raise ValueError('name already exists')
        self.short_names.append(name)
        if len(self.offsets):
            self.offsets.append(self.offsets[-1] + self.counts[-1])
        else:
            self.offsets.append(0)
        count = np.prod(dimensions)
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
    
    def pos(self, name, index):
        try:
            name_position = self.short_names.index(name)
        except ValueError:
            print ('name does not exist')
            raise
        offset = self.offsets[name_position]
        dimensions = self.dimensions_list[name_position]
        if np.isscalar(index) and np.isscalar(dimensions):
            return offset + index
        len_index = len(index)
        if len_index != len(dimensions):
            raise ValueError('index is the wrong shape')
        for i, num in enumerate(index):
            if num not in range(dimensions[i]):
                raise IndexError('index is out of range')
        if np.isscalar(index):
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


# TODO: this can probably be simplified - see pos method in variable class.
def generate_names_from_dimensions(name, dimensions):
    names = []
    if np.isscalar(dimensions):
        count = dimensions
        for n in range(dimensions):
            names.append("{}_{}".format(name, n))
    else:    
        periods = [1]
        ndims = len(dimensions)
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
    return names


# TODO: remove and replace with Variables.pos
def pos(variable, subscript):
    if variable == "x":
        if subscript[0] in range(N) and subscript[1] in range(P):
            return subscript[0] * P + subscript[1]
    elif variable == "y":
        if subscript[0] in range(M) and subscript[1] in range(P):
            return N * P + subscript[0] * P + subscript[1]
    elif variable == "z":
        if subscript in range(N):
            return N * P + M * P + subscript
    elif variable == "w":
        if (subscript[0] in range(N) and subscript[1] in range(N) and
            subscript[2] in range(P)):
            return (N * P + M * P + N + subscript[0] * P * N + subscript[1] * N
                    + subscript[2])
    elif variable == "v":
        if subscript[0] in range(N) and subscript[1] in range(N):
            return (N * P + M * P + N + N * N * P + subscript[0] * N 
                    + subscript[1])
    elif variable == "u":
        if subscript[0] in range(N) and subscript[1] in range(N):
            return (N * P + M * P + N + N * N * P + N * N + subscript[0] * N
                    + subscript[1])
    else:
        raise ValueError ("Invalid variable")
    raise ValueError ("Invalid subscript(s)")


def build_variables(parameters, buffers, vessels):
    # TODO: create function for initialising variables which feeds into tracker
    # of variables - possibly in a problem class, similarly create  a function
    # for initialising constraints, also in the same class
    # TODO: data should not be modulo cycle time, this should be handled here
    # somewhere and it should be possible to revert to original non-modulo data
        
    variables = Variables()
    
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

    # initialise indicator booleans u_(n,k): buffer k is used after buffer n
    # for all n in N for all k in N
    variables.add(name = "u",
                  dimensions = (N, N),
                  obj = [0] * N * N,
                  lb = [0] * N * N,
                  ub = [1] * N * N,
                  types = ["B"] * N * N)
    
    return variables

                       
def build_linear_constraints(parameters, buffers, vessels):
    c =  LinearConstraints()    
    counter = 0 # TODO: counter should live in the class, automatically update

    # TODO: define an add function in the class to tidy up all the items below.
    # each buffer can only be made in one slot
    for n in range(N):
        c.rhs.append(1)
        c.senses.append("E")
        c.names.append("one_slot_per_buffer_{}".format(n))
        for p in range(P):
            c.coefficients.append((counter, pos("x", (n, p)), 1))
        counter += 1
    
    # each slot contains a max of one vessel
    for p in range(P):
        c.rhs.append(1)
        c.senses.append("L")
        c.names.append("max_one_vessel_per_slot_{}".format(p))
        for m in range(M):
            c.coefficients.append((counter, pos("y", (m, p)), 1))
        counter += 1
    
    # vessel must be large enough
    for n in range(N):
        for p in range(P):
            c.rhs.append(0)
            c.senses.append("L")
            c.names.append("slot_{}_vessel_large_enough_for_buffer_{}".format(p, n))
            c.coefficients.append((counter, pos("x", (n, p)), buffers.volumes[n]))
            for m in range(M):
                c.coefficients.append((counter, pos("y", (m, p)), 0 - vessels.volumes[m]))
            counter += 1
    
    # vessel must be small enough
    # TODO: consider replacing with indicator - this may increase accuracy
    for n in range(N):
        for p in range(P):
            c.rhs.append(buffers.volumes[n] + max(vessels.volumes))
            c.senses.append("L")
            c.names.append("slot_{}_vessel_small_enough_for_buffer_{}".format(p, n))
            c.coefficients.append((counter, pos("x", (n, p)), max(vessels.volumes)))
            for m in range(M):
                c.coefficients.append((
                        counter, pos("y", (m, p)),
                        vessels.volumes[m] * parameters.minimum_fill_ratio))
            counter += 1
    
    # total duration in buffer hold must not exceed cycle time
    for n in range(N):
        c.rhs.append(
                parameters.cycle_time
                - parameters.hold_pre_duration
                - parameters.hold_post_duration
                - buffers.use_durations[n])
        c.senses.append("L")
        c.names.append(
                "buffer_{}_hold_vessel_duration_less_than_cycle_time".format(n))
        c.coefficients.append((counter, pos("z", n), 1))
        counter += 1
    
    # track if buffers n and k both made in slot p
    for n in range(N):
        for k in range(N):
            if k != n:
                for p in range(P):
                    c.rhs.append(1)
                    c.senses.append("L")
                    c.names.append(
                            "buffers_{}_and_{}_share_slot_{}".format(n, k, p))
                    c.coefficients.append((counter, pos("x", (n, p)), 1))
                    c.coefficients.append((counter, pos("x", (k, p)), 1))
                    c.coefficients.append((counter, pos("w", (n, k, p)), -1))
                    counter += 1
    
    # track if buffers n and k both made in same slot
    for n in range(N):
        for k in range(n + 1, N):
            c.rhs.append(0)
            c.senses.append("L")
            c.names.append("buffers_{}_and_{}_share_a_slot".format(n, k))
            c.coefficients.append((counter, pos("v", (n, k)), -1))
            for p in range(P):
                c.coefficients.append((counter, pos("w", (n, k, p)), 1))
            counter += 1
    
    # scheduling constraints - TODO: debug - doesn't seem to be working!!!
    prep_total_duration = (parameters.prep_pre_duration 
                           + parameters.transfer_duration
                           + parameters.prep_post_duration)
    for n in range(N):
        for k in range(n + 1, N):
            for i in (0, 1):
                c.rhs.append((2 * i - 1) * buffers.use_start_times[n]
                           + (1 - 2 * i) * buffers.use_start_times[k]
                           - prep_total_duration
                           + (i + 1) * 2 * parameters.cycle_time)
                c.senses.append("L")
                c.names.append("scheduling_buffer_{}_after_buffer_{}"
                             .format((k, n)[i], (n, k)[i]))
                c.coefficients.append((counter, pos("z", n), 2 * i - 1))
                c.coefficients.append((counter, pos("z", k), 1 - 2 * i))
                c.coefficients.append((counter, pos("u", (n, k)),
                                     2 * (2 * i - 1) * parameters.cycle_time))
                c.coefficients.append((counter, pos("v", (n, k)), 
                                     2 * parameters.cycle_time))
                counter += 1         
    
    return c


def add_variables(prob, variables):
    prob.variables.add(variables.obj, variables.lb, variables.ub, 
                       variables.types, variables.names)

def add_linear_constraints(prob, linear_constraints):
    prob.linear_constraints.add(rhs=linear_constraints.rhs,
                                senses=linear_constraints.senses,
                                names=linear_constraints.names)
    prob.linear_constraints.set_coefficients(linear_constraints.coefficients)


    
def add_indicator_constraints(prob):
    prep_total_duration = (parameters.prep_pre_duration 
                           + parameters.transfer_duration
                           + parameters.prep_post_duration)
       
    for n in range(N):
        for k in range(n + 1, N):
            indices = [pos("z", n),
                       pos("z", k),
                       pos("v", (n, k))]
            values = [[1, -1, 2 * parameters.cycle_time], 
                      [-1, 1, 2 * parameters.cycle_time]]
            for i in (0, 1):
                prob.indicator_constraints.add(
                        indvar=pos("v", (n, k)),
                        sense="L",
                        complemented=i,
                        lin_expr=[indices, values[i]],
                        name=("scheduling_buffers_{}_and_{}_"
                              "indicator_{}".format(n, k, i)),
                        rhs=(buffers.use_start_times[n] * (1 - 2 * i)
                             - buffers.use_start_times[k] * (1 - 2 * i)
                             - prep_total_duration 
                             + 2 * parameters.cycle_time))


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
def print_solution():
    values = prob.solution.get_values()
    selected_slots = []
    selected_volumes = []
    for m in range(M):
        for p in range(P):
            position = pos("y", (m,p))
            if abs(values[position] - 1.0) <= TOLERANCE:
                selected_slots.append(p)
                selected_volumes.append(vessels.volumes[m])
    lines = ["|slots    |"]
    for s in selected_slots: lines[0] += " {}\t\t|".format(s)
    lines.append("|volumes  |")
    for s in selected_volumes: lines[1] += " {}\t|".format(s)
    lines.append("\n|buffers  ")
    for n in range(N):
        line = "| {}\t  |".format(n)
        for p in selected_slots:
            position = pos("x", (n,p)) 
            if abs(values[position] - 1.0) <= TOLERANCE:
                line += " X\t\t|"
            else:
                line += "  \t\t|"
        lines.append(line)
    for l in lines: print(l)


# TODO: restructure the data using pandas


if __name__ == "__main__":
    parameters = Parameters()
    buffers = Buffers()
    vessels = Vessels()
    
    # TODO: where should these be???
    M = len(vessels.volumes)  # M = number of available vessel sizes
    N = len(buffers.volumes)  # N = number of buffers
    P = N  # P = number of slots for vessels (equal to number of buffers)
        
    linear_constraints = build_linear_constraints(parameters, buffers, vessels)
    variables = build_variables(parameters, buffers, vessels)
    
    prob = cplex.Cplex()
    prob.objective.set_sense(prob.objective.sense.minimize)
    add_variables(prob, variables)
    add_linear_constraints(prob, linear_constraints)
    add_indicator_constraints(prob)
    solve_model(prob)
