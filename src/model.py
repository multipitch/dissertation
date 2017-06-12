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
L = len(buffer_vol)
if len(time_of_first_use) != L or len(duration_of_use) != L:
    raise ValueError ("Buffer data incomplete")

# define vessel data
available_vessels = [1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000,
                     14000, 15000, 16000, 18000, 20000, 22000, 25000, 30000]
vessel_costs = [math.pow(i, 0.6) for i in available_vessels]
M  = len(available_vessels)

# define slots data
N = L

# initialise problem
prob = cplex.Cplex()    
prob.objective.set_sense(prob.objective.sense.minimize)

# initialise booleans x_l_n: buffer l made in slot n
for n in range(N):
    prob.variables.add(lb=[0.0]*L, ub=[1.0]*L, types=["B"]*L,
                       names=["x_{}_{}".format(l, n) for l in range(L)])

# initialise booleans y_m_n: vessel size m in slot n
for n in range(N):
    prob.variables.add(obj=vessel_costs, lb=[0.0]*M, ub=[1.0]*M, types=["B"]*M,
                       names=["y_{}_{}".format(m, n) for m in range(M)])

# initialise duration_hold_l: hold duration for all buffers l of L
prob.variables.add(lb=[duration_hold_min]*L, ub=[duration_hold_max]*L,
                   types=["C"]*L, 
                   names=["duration_hold_{}".format(l) for l in range(L)])

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
