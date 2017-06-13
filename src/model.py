import math
import cplex
from cplex.exceptions import CplexSolverError
#import numpy as np

#TODO: consider switching indices w and z???

# define general data TODO: repalce with data file
cycle_time = 4.0 * 24.0
prep_pre_duration = 8.0
prep_post_duration = 2.0
transfer_duration = 1.0
hold_pre_duration = 5.0
hold_post_duration = 2.0
hold_duration_min = 12.0
hold_duration_max = 36.0
minimum_fill_ratio = 0.3

# define buffer data TODO: replace with data file
buffer_volume = [12456.0, 18943.0, 2345.0, 6539.0, 21045.0, 7787.0, 10492.0]
buffer_use_start_time = [12.0, 32.0, 45.0, 27.0, 45.0, 45.0, 60.0]
buffer_use_duration = [2.0, 7.0, 20.0, 8.0, 20.0, 20.0, 12.0]
N = len(buffer_volume)  # N = number of buffers
if len(buffer_use_start_time) != N or len(buffer_use_duration) != N:
    raise ValueError ("Buffer data incomplete")

# define vessel data TODO: replace with data file
vessel_sizes = [
        1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 8000.0, 10000.0, 
        12000.0, 14000.0, 15000.0, 16000.0, 18000.0, 20000.0, 22000.0,
        25000.0, 30000.0]
vessel_costs = [math.pow(i, 0.6) for i in vessel_sizes]  # cost function
M  = len(vessel_sizes)  # M = number of available vessel sizes

# define slots data
P = N  # P = number of slots for vessels (equal to number of buffers)


# TODO: sanity check order of indices - have done some chopping and changing 
# and the results seem a bit strange as a result e.g. x_(n,p) vs x_(p,n) etc.

# initialise problem
prob = cplex.Cplex()    
prob.objective.set_sense(prob.objective.sense.minimize)

# initialise booleans x_(n,p): buffer n made in slot p for all n in N, p in P
prob.variables.add(
        obj=[0.0]*N*P, lb=[0.0]*N*P, ub=[1.0]*N*P, types=["B"]*N*P,
        names=["x_({},{})".format(n, p) for n in range(N) for p in range(P)])

# initialise booleans y_(m,p): vessel size m in slot p for all m in M, p in P
prob.variables.add(
        obj= [vessel_sizes[m] for m in range(M) for p in range(P)],
        lb=[0.0]*M*P, ub=[1.0]*M*P, types=["B"]*M*P,
        names=["y_({},{})".format(m, p) for m in range(M) for p in range(P)])


# initialise z_n: hold duration for buffer n for all n in N, incl bounds
prob.variables.add(
        obj=[0.0]*N, lb=[hold_duration_min]*N, ub=[hold_duration_max]*N, 
        types=["C"]*N, names=["z_{}".format(n) for n in range(N)])

# initialise indicator booleans w_(n,k,p): required for scheduling constraint, 
# for buffer n and buffer k and slot p for all n in N, k in N, p in P, k != n
nms = ["w_({},{},{})".format(n, k, p)
        for n in range(N) for k in range(N) for p in range(P)]
i = N*N*P
prob.variables.add(
        obj=[0.0]*i, lb=[0.0]*i, ub=[1.0]*i, types=["B"]*i, names=nms)
        

# define a function to keep track of the positions of the above variables
# TODO: rejigging indices might require this function to be broken - need to
# re-check
# TODO: do we need to add entry for w???
def pos(variable, subscript, M=M, N=N, P=P):
    if variable == "x":
        return subscript[1] * N + subscript[0]
    elif variable == "y":
        return N * P + subscript[1] * M + subscript[0]
    elif variable == "z":
        return N * P + M * P + subscript
    else:
        raise ValueError ("variable must be x or y or z")


# TODO: add constraints

# for debugging
names = prob.variables.get_names()

# for building constraints
right = []
types = []
coeff = []
cname = []
eqn = 0

# each buffer can only be made in one slot
eqn_offset = eqn
for n in range(N):
    right.append(1.0)
    types.append("E")
    cname.append("one_slot_per_buffer_{}".format(n))
    for p in range(P):
        coeff.append((eqn, pos("x", (n, p)), 1.0))
    eqn += 1    

# each slot contains a max of one vessel
eqn_offset = eqn
for p in range(P):
    right.append(1.0)
    types.append("L")
    cname.append("max_one_vessel_per_slot_{}".format(p))
    for m in range(M):
        coeff.append((eqn, pos("y", (m, p)), 1.0))
    eqn += 1    

# vessel must be large enough
eqn_offset = eqn
for n in range(N):
    for p in range(P):
        right.append(0.0)
        types.append("L")
        cname.append("vessel_large_enough_{}_{}".format(n, p))
        coeff.append((eqn, pos("x", (n, p)), buffer_volume[n]))
        for m in range(M):
            coeff.append((eqn, pos("y", (m, p)), 0.0 - vessel_sizes[m]))
        eqn += 1

"""
# vessel must be small enough
eqn_offset = eqn
for n in range(N):
    for p in range(P):
        right.append(0.0)
        types.append("L")
        cname.append("vessel_small_enough_{}_{}".format(n, p))
        coeff.append((eqn, pos("x", (n, p)), 0.0 - buffer_volume[n]))
        for m in range(M):
            coeff.append((
                    eqn, pos("y", (m, p)),
                    vessel_sizes[m] * minimum_fill_ratio))
        eqn += 1
"""

# total duration in buffer hold must not exceed cycle time
eqn_offset = eqn
for n in range(N):
    right.append(cycle_time - (hold_pre_duration + hold_post_duration +
            buffer_use_duration[n]))
    types.append("L")
    cname.append("hold_vessel_duration_less_than_cycle_time_{}".format(n))
    coeff.append((eqn, pos("z", n), 1.0))
    eqn += 1   

# apply above constraints to model
prob.linear_constraints.add(rhs=right, senses=types, names=cname)
prob.linear_constraints.set_coefficients(coeff)


# TODO: Add scheduling constraint (indicator function probably required)
# see /opt/ibm/ILOG/CPLEX_Studio1271/cplex/examples/src/python/etsp.py 


def equation_strings(coeff=coeff, types=types, right=right, cname=cname):
    equations = []
    for e in range(eqn):
        equation = ""
        first = True
        for c in coeff:
            if c[0] == e:
                if c[2] == 1:
                    equation += " + {}".format(names[c[1]])
                elif c[2] >= 0:
                    equation += " + {}{}".format(c[2], names[c[1]])
                else:
                    equation += " - {}{}".format(-c[2], names[c[1]])
                if first:
                    first = False
                    if c[2] >= 0:
                        equation = equation[3:]
        if types[e] == "L":
            equation += " <= "
        elif types[e] == "E":
            equation += " = "
        elif types[e] == "G":
            equation += " >= "
        elif types[e] == "R":
            raise TypeError("Shouldn't be using this type")
        else:
            raise ValueError ("Unknown type")
        equation += str(right[e])
        #equation = str(cname) + "\t:  " + equation
        equations.append(equation)
    return equations

equations = equation_strings()


# TODO: wrap the above in a function once it's working                      

# TODO: implement solver (based on included examples)

def solve_model():
    try:
        prob.solve()
    except CplexSolverError as e:
        print("Exception raised during solve: " + e)
    else:
        solution = prob.solution

        # solution.get_status() returns an integer code
        print("\nSolution status = ", solution.get_status(), ":", end=' ')
        # the following line prints the corresponding string
        print(solution.status[solution.get_status()])

        # Display solution.
        for i, j in enumerate(names): 
            print("{}\t:  {}".format(j, prob.solution.get_values(i)))


if __name__ == "__main__":
    solve_model()
