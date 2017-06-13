import math
import cplex
from cplex.exceptions import CplexError
import numpy as np

# define general data
cycle_time = 4 * 24
duration_prep_pre = 8
duration_prep_post = 2
duration_transfer = 1
duration_hold_pre = 5
duration_hold_post = 2
duration_hold_min = 12
duration_hold_max = 36
min_fill_ratio = 0.3

# define buffer data
buffer_vol = [12456, 18943, 2345, 6539, 21045, 7787, 10492]
time_of_first_use = [12, 32, 45, 27, 45, 45, 60]
duration_of_use = [2, 7, 20, 8, 20, 20, 12]
N = len(buffer_vol)
if len(time_of_first_use) != N or len(duration_of_use) != N:
    raise ValueError ("Buffer data incomplete")

# define vessel data
available_vessels = [1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000,
                     14000, 15000, 16000, 18000, 20000, 22000, 25000, 30000]
vessel_costs = [math.pow(i, 0.6) for i in available_vessels]
M  = len(available_vessels)

# define slots data
P = N

# initialise problem
prob = cplex.Cplex()    
prob.objective.set_sense(prob.objective.sense.minimize)

# initialise booleans x_n_p: buffer n made in slot p for all n in N, p in P
for p in range(P):
    prob.variables.add(
            lb=[0.0]*N, ub=[1.0]*N, types=["B"]*N,
            names=["x_{}_{}".format(n, p) for n in range(N)])

# initialise booleans y_m_p: vessel size m in slot p for all m in M, p in P
for p in range(P):
    prob.variables.add(
            obj=vessel_costs, lb=[0.0]*M, ub=[1.0]*M, types=["B"]*M,
            names=["y_{}_{}".format(m, p) for m in range(M)])

# initialise z_n: hold duration for buffer n for all n in N
prob.variables.add(
        lb=[duration_hold_min]*N, ub=[duration_hold_max]*N, types=["C"]*N,
        names=["z_{}".format(n) for n in range(N)])

# TODO: wrap the above in a function once it's working                      

# TODO: add constraints    
"""
prob.linear_constraints.add(lin_expr=rows, senses=my_sense,
                            rhs=my_rhs, names=my_rownames)
"""

# TODO: implement solver (based on included examples)
"""
try:
    model.solve()
except CplexSolverError as e:
    print("Exception raised during solve: " + e)
else:
    solution = model.solution

    # solution.get_status() returns an integer code
    print("Solution status = ", solution.get_status(), ":", end=' ')
    # the following line prints the corresponding string
    print(solution.status[solution.get_status()])

    # Display solution.
    print("Total cost = ", solution.get_objective_value())
"""

if __name__ == "__main__":
    pass
