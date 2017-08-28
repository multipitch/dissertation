import argparse
import configparser
import csv
import os
import random
import timeit

import numpy
import pulp

from plots import single_cycle_plot, timing_plot

# input parser
parser = argparse.ArgumentParser(prog="model.py", add_help=True,
                                 description="Run buffer vessel simulation.")
parser.add_argument("-P", "--path", default=".", type=str,
                    help="file path (default: <working directory>)")
parser.add_argument("-p", "--parameters", default="parameters.ini", type=str,
                    help="parameters filename (default: parameters.ini)")
parser.add_argument("-b", "--buffers", default="buffers.csv", type=str,
                    help="buffers filename (default: buffers.csv)")
parser.add_argument("-v", "--vessels", default="vessels.csv", type=str,
                    help="vessel filename (default: vessels.csv)")
parser.add_argument("-s", "--solver", default=None, type=str,
                    help="solver to be used (default: <PuLP default>)")
parser.add_argument("-n", "--no-secondary", action="store_true",
                    help="do not solve for secondary constraint")
parser.add_argument("-w", "--write", action="store_true",
                    help="write problem to file in .lp format")


# globals
ARGS = parser.parse_args()
PATH = os.path.abspath(os.path.expanduser(ARGS.path)) + "/"
PARAMETERS = PATH + ARGS.parameters
BUFFERS = PATH + ARGS.buffers
VESSELS = PATH + ARGS.vessels
WRITE = ARGS.write
NO_SECONDARY = ARGS.no_secondary
# the following solvers are installed on the test machine; solvers currently
# causing issues have been commented out
if ARGS.solver is None: # use pulp's default
    SOLVER = None 
elif ARGS.solver.upper() == "PULP_CBC_CMD":
    SOLVER = pulp.solvers.PULP_CBC_CMD(msg=1, threads=os.cpu_count())
elif ARGS.solver.upper() == "CPLEX":
    SOLVER = pulp.solvers.CPLEX()
elif ARGS.solver.upper() == "GLPK":
    SOLVER = pulp.solvers.GLPK()
#elif ARGS.solver.upper() == "XPRESS":
#    SOLVER = pulp.solvers.XPRESS()
#elif ARGS.solver.upper() in ["COIN", "COIN_CMD"]:
#    SOLVER = pulp.solvers.COIN_CMD()
else:
    raise ValueError("{} is an unsupported solver.".format(ARGS.solver))


class Parameters:
    
    def __init__(self, data=PARAMETERS):
        if type(data) == str:
            file_parser = configparser.SafeConfigParser()
            file_parser.read(data)
            parameters = file_parser.items("parameters")
            self.input_dict = dict([(p[0], float(p[1])) for p in parameters])
        elif type(data) == dict:
            self.input_dict = data
        else:
            raise TypeError("data should be a filename string or a dict")        
        self.cycle_time = self.input_dict["cycle_time"]
        self.prep_pre_duration = self.input_dict["prep_pre_duration"]
        self.prep_post_duration = self.input_dict["prep_post_duration"]
        self.transfer_duration = self.input_dict["transfer_duration"]
        self.hold_pre_duration = self.input_dict["hold_pre_duration"]
        self.hold_post_duration = self.input_dict["hold_post_duration"]
        self.hold_duration_min = self.input_dict["hold_duration_min"]
        self.hold_duration_max = self.input_dict["hold_duration_max"]
        self.minimum_fill_ratio = self.input_dict["minimum_fill_ratio"]
        self.maximum_prep_utilisation = (
                self.input_dict["maximum_prep_utilisation"])
        self.prep_total_duration = (self.prep_pre_duration
                                    + self.transfer_duration
                                    + self.prep_post_duration)


class Vessels:
    
    def __init__(self, data=VESSELS):
        if type(data) == str:
            self.input_dict = csv_columns_to_dict_of_lists(data)            
        elif type(data) == dict:
            self.input_dict = data
        else:
            raise TypeError("data should be a filename string or a dict")   
        self.names = self.input_dict["names"]
        self.volumes = self.input_dict["volumes"]
        self.costs = self.input_dict["costs"]
        self.max_volume = max(self.volumes)
        self.min_volume = min(self.volumes)
        self.count = len(self.names)


class Buffers:
    
    def __init__(self, cycle_time, data=BUFFERS):
        if type(data) == str:
            self.input_dict = csv_columns_to_dict_of_lists(data)        
        elif type(data) == dict:
            self.input_dict = data
        else:
            raise TypeError("data should be a filename string or a dict")  
        self.names = self.input_dict["names"]
        self.volumes = self.input_dict["volumes"]
        self.absolute_use_start_times = self.input_dict["use_start_times"]
        self.use_durations = self.input_dict["use_durations"]
        self.relative_use_start_times = []
        for t in self.absolute_use_start_times:
            self.relative_use_start_times.append(t % cycle_time)
        self.use_start_times = self.relative_use_start_times # shortened name
        self.count = len(self.names)
        self.prep_slots = None
        self.prep_vessels = None
        self.prep_start_times = None
        self.hold_start_times = None
        self.transfer_start_times = None
        self.prep_total_durations = None
        self.hold_total_durations = None
        self.transfer_durations = None
        self.p_to_m = None
    
    def get_results(self, parameters, results):
        ct = parameters.cycle_time
        n_to_p = dict(numpy.argwhere(results.x))
        self.p_to_m = dict(numpy.fliplr(numpy.argwhere(results.y)))
        p_to_m = self.p_to_m
        self.prep_slots = numpy.array([n_to_p[i] for i in range(self.count)])
        self.prep_vessels = numpy.array([p_to_m[i] for i in self.prep_slots])
        self.prep_start_times = (numpy.asarray(self.use_start_times)
                                 - results.z
                                 - parameters.transfer_duration
                                 - parameters.prep_pre_duration) % ct
        self.hold_start_times = (numpy.asarray(self.use_start_times)
                                 - results.z
                                 - parameters.transfer_duration
                                 - parameters.hold_pre_duration) % ct
        self.transfer_start_times = (self.prep_start_times
                                     + parameters.prep_pre_duration) % ct
        self.prep_total_durations = numpy.full(self.count,
                                                parameters.prep_total_duration)
        self.hold_total_durations = (parameters.hold_pre_duration
                                     + parameters.transfer_duration
                                     + results.z
                                     + numpy.asarray(self.use_durations)
                                     + parameters.hold_post_duration)
        self.transfer_durations = numpy.full(self.count,
                                             parameters.transfer_duration)


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
    for i in range(numpy.prod(dimensions)):
        varobject, position = next(variter)
        unflattened[position] = varobject
    return unflattened


def define_problem(parameters, buffers, vessels):
    M = vessels.count
    N = buffers.count
    P = N  # number of slots
    
    ct = parameters.cycle_time
    t_use = buffers.use_start_times
    t_prep = parameters.prep_total_duration
    
    # Problem
    problem = pulp.LpProblem("Buffer Preparation Vessel Selection",
                             pulp.LpMinimize)
    
    # Variables
    q = Variable("q", (N,), 0, 1, "Integer") # add ct to prep ref time
    r = Variable("r", (N,), 0, 1, "Integer") # add ct to free time lhs
    s = Variable("s", (N,), 0, 1, "Integer") # subtract ct from free time rhs
    u = Variable("u", (N,), 0, 1, "Integer") # select between free time cases
    v = Variable("v", (N, N), 0, 1, "Integer") # select between free time ranges
    w = Variable("w", (N, N, P), 0, 1, "Integer") # buffers n and k in slot p
    x = Variable("x", (N, P), 0, 1, "Integer") # buffer n in slot p
    y = Variable("y", (M, P), 0, 1, "Integer") # vessel m in slot p
    z = Variable("z", (N,), parameters.hold_duration_min,
                 parameters.hold_duration_max, "Continuous") # hold op duration
    
    # Keep track of variables in a class
    variables = Variables([q, r, s, u, v, w, x, y, z])
    
    # A buffer may only be prepared in one slot
    for n in range(N):
        problem += sum([x.o[n][p] for p in range(P)]) == 1
    
    # Each slot can contain at most one vessel
    for p in range(P):
        problem += sum([y.o[m][p] for m in range(M)]) <= 1
    
    # Vessels must be sized correctly for buffers
    mfr = parameters.minimum_fill_ratio
    vv = vessels.volumes
    for n in range(N):
        for p in range(P):
            problem += (buffers.volumes[n] * x.o[n][p]
                        <= sum([vv[m] * y.o[m][p] for m in range(M)]))
            problem += (vessels.max_volume * x.o[n][p]
                        + mfr * sum([vv[m]*y.o[m][p] for m in range(M)])
                        <= buffers.volumes[n] + vessels.max_volume)
    
    # In each slot, utilisation ratio must be below a maximum value
    mpu = parameters.maximum_prep_utilisation
    ptd = parameters.prep_total_duration
    for p in range(P):
        problem += sum([x.o[n][p] for n in range(N)]) * ptd <= ct * mpu
    
    # Hold procedure durations must be less than cycle time
    for n in range(N):
        rhs = (ct - parameters.hold_pre_duration - parameters.transfer_duration
               - parameters.hold_post_duration - buffers.use_durations[n])
        problem += z.o[n] <= rhs
    
    # For each pair of distinct buffers, for each slot, indicate if both
    # buffers are prepared in the given slot
    for n in range(N):
        for k in range(n + 1, N):
            for p in range(P):
                problem += x.o[n][p] + x.o[k][p] - w.o[n][k][p] <= 1
                problem += x.o[n][p] + x.o[k][p] >= 2 * w.o[n][k][p]
     
    
    # For each buffer, indicate if relative use start time minus hold time is
    # negative (defining q)
    for n in range(N):
        problem += ct * q.o[n] - z.o[n] >= - buffers.use_start_times[n]
        problem += ct * q.o[n] - z.o[n] <= - buffers.use_start_times[n] + ct
        
    # For each buffer, indicate if lower bound of free time window is greater
    # than the cycle time (defining r)
    for n in range(N):     
        problem += ct * r.o[n] - ct * q.o[n] + z.o[n] <= t_prep + t_use[n]
        problem += ct * r.o[n] - ct * q.o[n] + z.o[n] >= t_prep + t_use[n] - ct
    
    # For each buffer, indicate if upper bound of free time window is less
    # than zero (defining s)
    for n in range(N):
        problem += ct * q.o[n] + ct * s.o[n] - z.o[n] >= t_prep - t_use[n]
        problem += ct * q.o[n] + ct * s.o[n] - z.o[n] <= t_prep - t_use[n] + ct
    
    """
    # For each buffer, indicate if free time window crosses cycle boundary
    # (u = r xor s)
    for n in range(N):
        problem += r.o[n] + s.o[n] - u.o[n] >= 0
        problem += r.o[n] + s.o[n] + u.o[n] <= 2
        problem += r.o[n] - s.o[n] - u.o[n] <= 0
        problem += r.o[n] - s.o[n] + u.o[n] >= 0
    """
    
    # For each buffer, indicate if free time window crosses cycle boundary
    # (u = r or s)
    for n in range(N):
        problem += r.o[n] + s.o[n] - u.o[n] >= 0
        problem += r.o[n] + s.o[n] - 2 * u.o[n] <= 0
    
    # Each prep vessel can only do one thing at a time
    big_M = 2 * ct
    for n in range(N):
        for k in range(n + 1, N):
            problem += (ct * q.o[k] - ct * q.o[n] - z.o[k] + z.o[n]
                        + big_M * u.o[n] + big_M * v.o[n][k]
                        - big_M * sum([w.o[n][k][p] for p in range(P)])
                        >= t_use[n] - t_use[k] + t_prep - big_M)
            problem += (ct * q.o[k] - ct * q.o[n] - z.o[k] + z.o[n]
                        - big_M * u.o[n] + big_M * v.o[n][k]
                        + big_M * sum([w.o[n][k][p] for p in range(P)])
                        <= t_use[n] - t_use[k] - t_prep + 2 * big_M)
            problem += (ct * q.o[k] - ct * q.o[n] - z.o[k] + z.o[n]
                        - big_M * u.o[n] + ct * s.o[n]
                        - big_M * sum([w.o[n][k][p] for p in range(P)])
                        >= t_use[n] - t_use[k] + t_prep - 2 * big_M)
            problem += (ct * q.o[k] - ct * q.o[n] - z.o[k] + z.o[n]
                        + big_M * u.o[n] - ct * r.o[n]
                        + big_M * sum([w.o[n][k][p] for p in range(P)])
                        <= t_use[n] - t_use[k] - t_prep + 2 * big_M)
    
    """
    # Hack to avoid edge cases at upper cycle time boundary
    # TODO: May not be required - test robustness with this disabled
    for n in range(N):
        problem += ct * q.o[n] - z.o[n] <= - t_use[n] + (1 - DEADBAND) * ct
        problem += (ct * q.o[n] + ct * r.o[n] - z.o[n]
                    <= t_prep - t_use[n] + (1 - DEADBAND) * ct)
        problem += (ct * q.o[n] - ct * s.o[n] - z.o[n] 
                    <= - t_prep - t_use[n] + (1 - DEADBAND) * ct)
    """
      
    return problem, variables
    

def initial_objective(problem, variables, buffers, vessels):
    M = vessels.count
    N = buffers.count
    P = N  # number of slots
    
    # Minimise total vessel cost
    problem += sum([vessels.costs[m] * variables.y.o[m][p] 
                    for m in range(M) for p in range(P)])


def secondary_objective(problem, variables, buffers, vessels,
                        initial_objective_value):
    M = vessels.count
    N = buffers.count
    P = N  # number of slots
    
    # Constrain total vessel cost to be equal to initial objective
    vc = vessels.costs
    ys = variables.y.o
    problem += (sum([vc[m] * ys[m][p] for m in range(M) for p in range(P)])
                == initial_objective_value)
    
    # Minimise sum of hold times
    problem += sum([variables.z.o[n] for n in range(N)])


def generate_random_model(N, min_duration_ratio=0.2, max_duration_ratio=0.9,
                          to_file=False):
    parameters = Parameters()
    vessels = Vessels()
    buffer_dict = {"names": [], "volumes": [], "use_start_times": [],
                   "use_durations": []}
                   
    ct = parameters.cycle_time
    max_vol = vessels.max_volume
    min_vol = vessels.min_volume * parameters.minimum_fill_ratio
    min_duration = min_duration_ratio * ct
    max_feasible_duration = ct - (parameters.hold_pre_duration
                                  + parameters.transfer_duration
                                  + parameters.hold_duration_min
                                  + parameters.hold_post_duration)
    max_duration = max_feasible_duration * max_duration_ratio    
    
    for n in range(N):
        buffer_dict["names"].append("Buffer #{}".format(n + 1))  
        volume = round(random.uniform(min_vol, max_vol), 2)
        buffer_dict["volumes"].append(volume)
        use_start_time = round(random.uniform(0, ct), 2)
        buffer_dict["use_start_times"].append(use_start_time)
        use_duration = round(random.uniform(min_duration, max_duration), 2)
        buffer_dict["use_durations"].append(use_duration)
       
    if to_file:
        with open("buffers.csv", "w") as f:
            w = csv.writer(f, delimiter=",", quoting=csv.QUOTE_NONNUMERIC,
                           quotechar="\"",)
            keys = ["names", "volumes", "use_start_times", "use_durations"]
            w.writerow(keys)
            w.writerows(zip(*[buffer_dict[key] for key in keys]))
            return

    buffers = Buffers(parameters.cycle_time, buffer_dict)
    return parameters, vessels, buffers

    
def run_primary(solver=SOLVER, plot=True, write=False):
    parameters = Parameters()
    buffers = Buffers(parameters.cycle_time)
    vessels = Vessels() 
    problem, variables = define_problem(parameters, buffers, vessels)    
    # Optimise for initial objective: minimise cost
    initial_objective(problem, variables, buffers, vessels)
    if write:
        problem.writeLP("primary.lp")
    problem.solve(solver)
    print("Status:", pulp.LpStatus[problem.status])
    initial_objective_value = pulp.value(problem.objective)
    results = Results(variables)
    buffers.get_results(parameters, results)
    if plot:
        single_cycle_plot(parameters, buffers, vessels, (PATH + "plot1.pdf"))
    return (parameters, buffers, vessels, problem, variables,
            initial_objective_value)


def run_secondary(parameters, buffers, vessels, problem, variables,
                  initial_objective_value, solver=SOLVER, plot=True,
                  write=False):
    # Optimise for secondary objective: minimise hold times
    secondary_objective(problem, variables, buffers, vessels,
                        initial_objective_value)
    if write:
        problem.writeLP("secondary.lp")
    problem.solve(solver)
    print("Status:", pulp.LpStatus[problem.status])
    secondary_objective_value = pulp.value(problem.objective)
    results = Results(variables)
    buffers.get_results(parameters, results)
    if plot:
        single_cycle_plot(parameters, buffers, vessels, (PATH + "plot2.pdf"))
    return (parameters, buffers, vessels, problem, variables,
            secondary_objective_value)


def standard_run(no_secondary=NO_SECONDARY, write=WRITE, solver=SOLVER):
        primary = run_primary(solver=solver, plot=True, write=write)
        if not no_secondary:
            secondary = run_secondary(*primary, solver=solver, plot=True, 
                                      write=write)


def run_random_models(sizes, count=100, to_file=True, verbose=True, ret=True):
    #random.seed(3876401295) # for repeatability
    durations = numpy.empty((len(sizes), count))
    
    if to_file:
        with open("durations.csv", "w") as f:
            writer = csv.writer(f)
            writer.writerow(sizes)

    for n, N in enumerate(sizes):
        for c in range(count):
            parameters, vessels, buffers = generate_random_model(N)
            problem, variables = define_problem(parameters, buffers, vessels)
            initial_objective(problem, variables, buffers, vessels)
            start = timeit.default_timer()
            try:
                problem.solve(SOLVER)
                end = timeit.default_timer()
                if problem.status != 1:
                    durations[n][c] = numpy.inf
                    print("\n\nWarning: Optimium not found for run {}, size {}"
                          .format(c, N))
                else:
                    duration = end - start
                    durations[n][c] = duration
                if verbose:
                    print("\n\nCompleted run {}, size {}".format(c, N))
                print("Solver time: {0:.32f} s\n\n".format(duration))
            except:
                durations[n][c] = numpy.inf
                print("\n\nWarning: Optimium not found for run {}, size {}"
                      .format(c, N))
              
        if to_file:
            with open("durations.csv", "ab") as f:
                numpy.savetxt(f, durations, delimiter=",")

    if ret:
        return durations


def many_random(N):
    for n in range(N):
        generate_random_model(12, min_duration_ratio=0.2,
                              max_duration_ratio=0.9, to_file=False)
        primary = run_primary(plot=False)
       

def one_random(seed=None):
    if seed:
        random.seed(seed)
    generate_random_model(12, min_duration_ratio=0.2,
                          max_duration_ratio=0.9, to_file=True)
    standard_run()

if __name__ == "__main__":
    # TODO: error handling for failed optimisations
    
    #standard_run()
   
    #many_random(100)
    
    # the line below was used to generate example data used in ch3/ch4
    #one_random(123456)    
    
    #sizes = [2, 4, 6, 8, 10, 12]
    #durations = run_random_models(sizes, 100)
    #avg_durations = dict(zip(sizes, durations.mean(axis=1)))
    #timing_plot(sizes, durations)
    
    pass
