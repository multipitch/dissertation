import csv
import math

import cplex
import pandas as pd
import numpy as np

from configparser import SafeConfigParser

from cplex.exceptions import CplexSolverError


# define global parameters
TOLERANCE = 1E-6

# TODO: create nice interface so that constraints can be switched off/on for a
# run


class Parameters():
    
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


class DataClass():

    def __init__(self, filename):
        self.filename = filename        
        self.__dict__ = csv_columns_to_dict_of_lists(self.filename)

    def __repr__(self):
        return str(self.__dict__)
        

class Vessels(DataClass):
    
    def __init__(self, filename="vessels.csv"):
        self.filename = filename
        self.__dict__ = csv_columns_to_dict_of_lists(self.filename)


class Buffers(DataClass):
    
    def __init__(self, filename="buffers.csv"):
        self.filename = filename
        self.__dict__ = csv_columns_to_dict_of_lists(self.filename)


def csv_columns_to_dict_of_lists(filename):    
    with open(filename) as f:
        reader = csv.reader(f, skipinitialspace=True, delimiter=",",
                            quoting=csv.QUOTE_NONNUMERIC)
        data_list = list(reader)
        lines = len(data_list)
        data_dict = {}        
        for i, key in enumerate(data_list[0]):
            data_dict[key] = [data_list[j][i] for j in range(1, lines)]
        return data_dict
                

def build_model(params, buffers, vessels):
    # TODO: create function for initialising variables which feeds into tracker of
    # variables - possibly in a problem class, similarly create  a function for 
    # initialising constraints, also in the same class

    prob = cplex.Cplex()
    prob.objective.set_sense(prob.objective.sense.minimize)

    buffer_volumes = buffers.volumes
    buffer_use_start_times = buffers.use_start_times
    buffer_use_durations = buffers.use_durations
    vessel_volumes = vessels.volumes
    vessel_costs = vessels.costs

    
    # derive some values from the data
    M  = len(vessel_volumes)  # M = number of available vessel sizes
    N = len(buffer_volumes)  # N = number of buffers
    P = N  # P = number of slots for vessels (equal to number of buffers)
    prep_total_duration = (
            params.prep_pre_duration 
            + params.transfer_duration
            + params.prep_post_duration)
    
    # check buffer data
    if len(buffer_use_start_times) != N or len(buffer_use_durations) != N:
        raise ValueError ("Buffer data incomplete")
    for n in buffer_use_start_times:
        if n >= params.cycle_time:
            raise ValueError("buffer_use_start_times must be modulo cycle time")
    for n in buffer_volumes:
        if n > max(vessel_volumes):
            raise ValueError("buffer volume(s) greater than max vessel volume")
            
    # initialise booleans x_(n,p): buffer n made in slot p for all n in N, p in P
    prob.variables.add(
            obj=[0]*N*P, lb=[0]*N*P, ub=[1]*N*P, types=["B"]*N*P,
            names=["x_({},{})".format(n, p) for n in range(N) for p in range(P)])

    # initialise booleans y_(m,p): vessel size m in slot p for all m in M, p in P
    prob.variables.add(
            obj= [vessel_volumes[m] for m in range(M) for p in range(P)],
            lb=[0]*M*P, ub=[1]*M*P, types=["B"]*M*P,
            names=["y_({},{})".format(m, p) for m in range(M) for p in range(P)])


    # initialise z_n: hold duration for buffer n for all n in N, incl the following
    # bounds: params.hold_duration_min <= z_n <= params.hold_duration_max
    prob.variables.add(
            obj=[0]*N, lb=[params.hold_duration_min]*N, ub=[params.hold_duration_max]*N,
            types=["C"]*N, names=["z_{}".format(n) for n in range(N)])

    # initialise indicator booleans w_(n,k,p): buffers n and k both made in slot p
    inames = [
            "w_({},{},{})".format(n, k, p) for n in range(N) for k in range(N)
            for p in range(P)]
    i = N * N * P
    prob.variables.add(
            obj=[0]*i, lb=[0]*i, ub=[1]*i, types=["B"]*i, names=inames)

    # initialise indicator booleans v_(n,k): buffers n and k share a slot
    inames = ["v_({},{})".format(n, k) for n in range(N) for k in range(N)]
    i = N * N
    prob.variables.add(
            obj=[0]*i, lb=[0]*i, ub=[1]*i, types=["B"]*i, names=inames)


    # define a function to keep track of the positions of the above variables
    def pos(variable, subscript, M=M, N=N, P=P):
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
            if (
                    subscript[0] in range(N) and subscript[1] in range(N) and
                    subscript[2] in range(P)):
                return (
                        N * P + M * P + N + subscript[0] * P * N +
                        subscript[1] * N + subscript[2])
        elif variable == "v":
            if subscript[0] in range(M) and subscript[1] in range(P):
                return (
                        N * P + M * P + N + N * N * P + subscript[0] * N +
                        subscript[1])
        else:
            raise ValueError ("Invalid variable")
        raise ValueError ("Invalid subscript(s)")


    # for debugging
    vnames = prob.variables.get_names()

    # for building constraints
    rights = []
    types = []
    coeffs = []
    cnames = []
    rangevals = []
    eqn = 0

    # each buffer can only be made in one slot
    for n in range(N):
        rights.append(1)
        types.append("E")
        cnames.append("one_slot_per_buffer_{}".format(n))
        rangevals.append(0)
        for p in range(P):
            coeffs.append((eqn, pos("x", (n, p)), 1))
        eqn += 1

    # each slot contains a max of one vessel
    for p in range(P):
        rights.append(1)
        types.append("L")
        cnames.append("max_one_vessel_per_slot_{}".format(p))
        rangevals.append(0)
        for m in range(M):
            coeffs.append((eqn, pos("y", (m, p)), 1))
        eqn += 1

    # vessel must be large enough
    for n in range(N):
        for p in range(P):
            rights.append(0)
            types.append("L")
            cnames.append("slot_{}_vessel_large_enough_for_buffer_{}".format(p, n))
            rangevals.append(0)
            coeffs.append((eqn, pos("x", (n, p)), buffer_volumes[n]))
            for m in range(M):
                coeffs.append((eqn, pos("y", (m, p)), 0 - vessel_volumes[m]))
            eqn += 1

    # vessel must be small enough
    # TODO: consider replacing with indicator - this may increase accuracy
    for n in range(N):
        for p in range(P):
            rights.append(buffer_volumes[n] + max(vessel_volumes))
            types.append("L")
            cnames.append("slot_{}_vessel_small_enough_for_buffer_{}".format(p, n))
            rangevals.append(0)
            coeffs.append((eqn, pos("x", (n, p)), max(vessel_volumes)))
            for m in range(M):
                coeffs.append((
                        eqn, pos("y", (m, p)),
                        vessel_volumes[m] * params.minimum_fill_ratio))
            eqn += 1

    # total duration in buffer hold must not exceed cycle time
    for n in range(N):
        rights.append(
                params.cycle_time
                - params.hold_pre_duration
                - params.hold_post_duration
                - buffer_use_durations[n])
        types.append("L")
        cnames.append(
                "buffer_{}_hold_vessel_duration_less_than_cycle_time".format(n))
        rangevals.append(0)
        coeffs.append((eqn, pos("z", n), 1))
        eqn += 1

    # track if buffers n and k both made in slot p
    for n in range(N):
        for k in range(N):
            if k != n:
                for p in range(P):
                    rights.append(1)
                    types.append("L")
                    cnames.append(
                            "buffers_{}_and_{}_share_slot_{}".format(n, k, p))
                    rangevals.append(0)
                    coeffs.append((eqn, pos("x", (n, p)), 1))
                    coeffs.append((eqn, pos("x", (k, p)), 1))
                    coeffs.append((eqn, pos("w", (n, k, p)), -1))
                    eqn += 1

    # track if buffers n and k both made in same slot
    for n in range(N):
        for k in range(n + 1, N):
            rights.append(0)
            types.append("L")
            cnames.append("buffers_{}_and_{}_share_a_slot".format(n, k))
            rangevals.append(0)
            coeffs.append((eqn, pos("v", (n, k)), -1))
            for p in range(P):
                coeffs.append((eqn, pos("w", (n, k, p)), 1))
            eqn += 1

    # scheduling constraints (indicator for absolute time differences, with a big-M
    # approach for activating only if buffers n and k made in same slot
    for n in range(N):
        for k in range(n + 1, N):
            indices = [pos("z", n), pos("z", k), pos("v", (n, k))]
            values = [
                    [1, -1, 2 * params.cycle_time], 
                    [-1, 1, 2 * params.cycle_time]]
            for i in (0, 1):
                prob.indicator_constraints.add(
                        indvar=pos("v", (n, k)),
                        sense="L",
                        complemented=i,
                        lin_expr=[indices, values[i]],
                        name=(
                                "scheduling_buffers_{}_and_{}_"
                                "indicator_{}".format(n, k, p, i)),
                        rhs=(
                                buffer_use_start_times[n] * (1 - 2 * i) -
                                buffer_use_start_times[k] * (1 - 2 * i) -
                                prep_total_duration + 2 * params.cycle_time))


    #TODO: include strings for indicator constraints in function below
    #TODO: consider removing or putting in a seperate file
    def constraint_strings(
            count=eqn, coeffs=coeffs, types=types, rights=rights, cnames=cnames):
        constraint_strings = []
        for constraint in range(count):
            constraint_string = ""
            first = True
            for c in coeffs:
                if c[0] == constraint:
                    if c[2] == 1:
                        constraint_string += " + {}".format(vnames[c[1]])
                    elif c[2] >= 0:
                        constraint_string += " + {}{}".format(c[2], vnames[c[1]])
                    else:
                        constraint_string += " - {}{}".format(-c[2], vnames[c[1]])
                    if first:
                        first = False
                        if c[2] >= 0:
                            constraint_string = constraint_string[3:]
                        else:
                            constraint_string = "-" + constraint_string[3:]
            if types[constraint] == "L":
                constraint_string += " <= "
            elif types[constraint] == "E":
                constraint_string += " = "
            elif types[constraint] == "G":
                constraint_string += " >= "
            elif types[constraint] == "R":
                constraint_string += " TBC "
            else:
                raise TypeError ("Unknown type")
            constraint_string += str(rights[constraint])
            constraint_strings.append(constraint_string)
        return constraint_strings


    constraints = constraint_strings()


    #TODO: consider removing or putting in a seperate file
    def print_constraints():
        print(
                "w_(n, k, p):\tboolean  indicator variable for ensuring buffers n "
                "and k do not \n\t\t\t clash in slot p (boolean) for n in N, for"
                " k in N, for\n\t\t\t p in P")
        print("x_(n,p):\tboolean  buffer n made in slot p, for n in N, for p in P")
        print("y_(m,p):\tboolean  vessel size m in slot p, for m in M, for p in P")
        print("z_n:\t\treal\t buffer n hold duration, for n in N")
        for i, j in enumerate(constraints):
            print("{}\t:  {}".format(cnames[i], j))


    # apply above constraints to model
    prob.linear_constraints.add(
            rhs=rights, senses=types, names=cnames, range_values=rangevals)
    prob.linear_constraints.set_coefficients(coeffs)


    # TODO: wrap all the above in functions once it's working
    return prob

def solve_model(prob):
    try:
        prob.solve()
    except CplexSolverError as e:
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
                selected_volumes.append(vessel_volumes[m])
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
    params = Parameters()
    buffers = Buffers()
    vessels = Vessels()
    print(params)
    print(buffers)
    print(vessels)
    #buffers = pd.read_csv("buffers.csv")
    #vessels = pd.read_csv("vessels.csv")
    #prob = initialise_problem()
    prob = build_model(params, buffers, vessels)
    solve_model(prob)
