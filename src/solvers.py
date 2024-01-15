from utils import (
    PTN,
    EAN,
    pesp_part_to_betas,
    is_positive_definite,
    pesp_weights_to_betas,
    )
import pandas as pd
import gurobipy as gp
import numpy as np
import logging
import csv
import time as time_module


class ColppSolver_beta:
    '''
    Class for solving the cost oriented line planning problem with binary
        frequency decision variables. If input betas are none, then it is the regular COLPP solver.
    '''

    def __init__(
            self,
            ptn: PTN,
            linepool: dict,
            beta=None,
            test_lineplan=None,
            run_solve=True,
            TIMELIMIT=None,
            NODEFILESTART=None) -> None:

        logging.info('INIT COLPP')
        # input variables
        self.ptn = ptn
        self.linepool = linepool
        self.beta = beta
        self.test_lineplan = test_lineplan

        # output variables
        self.lineplan = None
        self.obj_val = None

        # optimization params
        self.TIMELIMIT = TIMELIMIT
        self.NODEFILESTART = NODEFILESTART

        # Solve
        if run_solve:
            self.solve()

    def solve(self) -> None:
        '''Solves the line planning problem'''

        # notation
        # V = self.ptn.V # not acessed
        E = self.ptn.E

        logging.info('Initialize model')
        model = gp.Model('COLPP')

        logging.info('Initialize decision variables')
        f = {l: model.addVar(vtype=gp.GRB.BINARY, name=f'freq_{l}')
             for l in self.linepool}

        # optional: test a line plan by fixing variables
        if self.test_lineplan is not None:
            logging.info('Fixing test line plan variables')
            self.test_lineplan = {
                key.replace('freq_', ''): value for (key, value) in self.test_lineplan.items()}
            for line in self.test_lineplan:
                if line in self.linepool:
                    f[line].start = self.test_lineplan[line]
                    f[line].lb = self.test_lineplan[line]
                    f[line].ub = self.test_lineplan[line]

        # objective function
        logging.info('Defining objective function')
        cost = gp.quicksum(
            f[line] * self.linepool[line]['cost'] for line in self.linepool)

        if self.beta is not None:
            label_to_index = {key: index for index, key in enumerate(self.linepool)}
            cost += gp.quicksum(
                self.beta[label_to_index[line1], label_to_index[line2]] * f[line1] * f[line2]
                for line1 in self.linepool
                for line2 in self.linepool
            )

        model.setObjective(cost, gp.GRB.MINIMIZE)

        logging.info('Adding constraints')
        # dictionary with edges as keys and a list of lines which operate
        #   the edge
        edges_of_line = {edge: [] for edge in E}
        for line in self.linepool:
            stations = self.linepool[line]['stations']
            edges_between_stations = zip(
                stations, stations[1:]
            )
            for edge in edges_between_stations:
                if edge in edges_of_line:  # key error for some reason
                    edges_of_line[edge].append(f[line])
                else:
                    print(f'{edge} not found in PTN.E?!')

        # Now add for each edge the constraint the minimal and maximal frequency
        for edge in edges_of_line:
            model.addConstr(
                gp.quicksum(edges_of_line[edge]) >= E[edge]['fmin'],
                name=f'fmin_{edge}'
            )
            model.addConstr(
                gp.quicksum(edges_of_line[edge]) <= E[edge]['fmax'],
                name=f'fmax_{edge}'
            )

        # save storage space
        del edges_of_line

        # setting optimization parameters
        if self.TIMELIMIT is not None:
            model.setParam('TimeLimit', self.TIMELIMIT)

        if self.NODEFILESTART is not None:
            model.setParam('NodefileStart', self.NODEFILESTART)

        # optimize the model
        logging.info('Solving the model')
        model.optimize()

        # read out optimal solution
        if model.status == gp.GRB.INFEASIBLE:
            logging.warning('Model is infeasible!')
            logging.info('Calculating smallest set of infeasible constraints')
            model.computeIIS()
            model.write(f'model_{model.modelName}.ilp')

        elif model.status == gp.GRB.OPTIMAL or model.status == gp.GRB.TIME_LIMIT:
            logging.info(f'Model terminated with the status {model.status}!')
            var_values = {}
            all_vars = model.getVars()
            var_values = {var.varName: var.x for var in all_vars}
            self.lineplan = var_values
            self.obj_val = model.objVal

        return

    def print_optimal_solution(self) -> None:

        print(f'Objective value: {self.obj_val}')
        print('Optimal variable values:')
        for name, val in self.lineplan.items():
            print(f'{name} = {val}')

        return


class PespSolver:
    '''Class for solving the periodic event scheduling problem.'''
    def __init__(
            self,
            ean: EAN,
            T=60,
            test_schedule=None,
            run_solve=True,
            TIMELIMIT=None,
            NODEFILESTART=None) -> None:

        # input Var
        logging.info('INIT PESP')
        self.ean = ean
        self.T = T
        self.test_schedule = test_schedule

        # output var
        self.schedule = None
        self.obj_val = None
        self.z_modulo_param = None

        # optimization params
        self.TIMELIMIT = TIMELIMIT
        self.NODEFILESTART = NODEFILESTART

        # solve
        if run_solve:
            self.solve()

    def solve(self) -> None:
        '''Solves the periodic event scheduling problem for a given event
        activity network with period time T'''

        # notation
        Events = self.ean.get_unique_events()
        Activities = self.ean.Activities

        # initialize model
        logging.info('Initialize model')
        model = gp.Model('PESP')

        # initialize decision variables
        logging.info('Initialize decision variables')
        pi = {}
        z = {}

        # pi - decisionvariable
        for event in Events:
            pi[event] = model.addVar(
                lb=0,
                ub=self.T-1,
                vtype=gp.GRB.INTEGER,
                name=f'pi_{event}')
            
        # z - decisionvariable
        for activity in Activities:
            z[activity] = model.addVar(
                vtype=gp.GRB.INTEGER, name=f'z_{activity}')

        # optional: test a schedule by fixing variables
        if self.test_schedule is not None:
            logging.info('Fixing test schedule variables')
            for event in self.test_schedule:
                if event in Events:
                    pi[event].start = self.test_schedule[event]
                    pi[event].lb = self.test_schedule[event]
                    pi[event].ub = self.test_schedule[event]

        # objective function
        logging.info('Defining objective function')
        duration = 0
        for activity, properties in self.ean.Activities.items():
            weight = properties['weight']
            i, j = activity[0], activity[1]
            duration += weight * (pi[j] - pi[i] + self.T * z[activity])

        model.setObjective(duration, gp.GRB.MINIMIZE)

        # constraints
        logging.info('Adding constraints')
        for activity, properties in self.ean.Activities.items():
            lb, ub = properties['lb'], properties['ub']
            i, j = activity[0], activity[1]
            model.addConstr(
                pi[j] - pi[i] + self.T * z[activity] >= lb,
                name=f'lb_{(i, j)}')
            model.addConstr(
                pi[j] - pi[i] + self.T * z[activity] <= ub,
                name=f'ub_{(i, j)}')
            
        # setting optimization parameters
        if self.TIMELIMIT is not None:
            model.setParam('TimeLimit', self.TIMELIMIT)

        if self.NODEFILESTART is not None:
            model.setParam('NodefileStart', self.NODEFILESTART)

        # optimize the model
        # model.write(f'pre_optimize.lp')
        logging.info('Solving the model')
        model.optimize()

        # read out optimal solution
        if model.status == gp.GRB.INFEASIBLE:
            logging.warning('Model is infeasible!')
            logging.info('Calculating smallest set of infeasible constraints')
            model.computeIIS()
            model.write(f'model_{model.modelName}.ilp')

        elif model.status == gp.GRB.OPTIMAL or model.status == gp.GRB.TIME_LIMIT:
            logging.info(f'Model terminated with the status {model.status}!')
            self.obj_val = model.objVal
            self.schedule = {event: pi[event].x for event in Events}
            self.z_modulo_param = {
                activity: z[activity].x for activity in Activities}

        else:
            print('SHOULD NOT HAPPEN No optimal solution found.')

        return
    
    def print_sorted_optimal_schedule(self) -> None:
        sorted_schedule = sorted(
            self.schedule.items(),
            key=lambda x: (x[0][0], x[0][2], x[0][1]))
        print(f'Objective value: {self.obj_val}')
        print('Optimal variable values:')
        for event, schedule_time in sorted_schedule:
            print(f'{event}: {schedule_time}')
        return
    
    def print_optimal_schedule(self) -> None:

        print(f'Objective value: {self.obj_val}')
        print('Optimal variable values:')
        for name, val in self.schedule.items():
            print(f'{name} = {val}')

        return
    

class LintimSolver:
    '''
    Solves the integrated line planning and timetabling problem for a given
    public transport network, event-activity network and a given time period T.
    '''

    def __init__(
            self,
            ptn: PTN,
            ean: EAN,
            linepool: dict,
            T: int,
            test_lineplan=None,
            test_schedule=None,
            start_pi=None,
            start_z=None,
            start_y=None,
            start_x=None,
            run_solve=True,
            TIMELIMIT=None,
            SAVE_CALLBACK=False) -> None:

        logging.info('INIT LinTim')
        # in- and output vars: COLPP
        self.ptn = ptn
        self.linepool = linepool
        self.test_lineplan = test_lineplan
        self.lineplan = None

        # in- and output vars: PESP
        self.ean = ean
        self.T = T
        self.test_schedule = test_schedule
        self.schedule = None
        self.obj_val = None
        self.z_modulo_param = None

        # additional variables
        self.M = self.T + 1
        self.y_connection_param = None
        self.x_connection_param = None
        self.obj_val = None
        self.SAVE_CALLBACK = SAVE_CALLBACK
        self.start_pi = start_pi
        self.start_z = start_z
        self.start_y = start_y
        self.start_x = start_x

        # optimization params
        self.TIMELIMIT = TIMELIMIT

        # solve
        if run_solve:
            self.solve()

    def solve(self) -> None:
            '''Solves the integrated line planning and timetabling problem for a
                period time T'''
            
            # for logging time - objVal
            def my_callback(model, where):
                if where == gp.GRB.Callback.MIP:
                    # Access and store the current best objective value and runtime
                    obj_val = model.cbGet(gp.GRB.Callback.MIP_OBJBST)
                    current_time = time_module.time() - start_time
                    # Store in a list
                    callback_data.append((current_time, obj_val))


            # notation
            # V = self.ptn.V  # not acessed
            E = self.ptn.E
            Events = self.ean.get_unique_events()
            Activities = self.ean.Activities

            # initialize model
            logging.info('Initialize model')
            model = gp.Model('Integrated_Model')
                
            # initialize decision variables
            logging.info('Initialize decision variables')

            # COLPP decision variables
            f = {l: model.addVar(vtype=gp.GRB.BINARY, name=f'freq_{l}') 
                for l in self.linepool}
            
            # PESP decision variables
            pi = {}
            z = {}

            for event in Events:
                pi[event] = model.addVar(
                    lb=0,
                    ub=self.T-1,
                    vtype=gp.GRB.INTEGER,
                    name=f'pi_{event}')
                
            for activity in Activities:
                z[activity] = model.addVar(
                    vtype=gp.GRB.INTEGER, name=f'z_{activity}')

            # Connection variables
            y = {}
            x = {}

            for activity in Activities:
                y[activity] = model.addVar(
                    vtype=gp.GRB.BINARY, name=f'y_{activity}')
                x[activity] = model.addVar(
                    lb=0, vtype=gp.GRB.INTEGER, name=f'x_{activity}')

            # optional: test a line plan by fixing variables
            if self.test_lineplan is not None:
                logging.info('Fixing test line plan variables')
                self.test_lineplan = {
                    key.replace('freq_', ''): value for (key, value) in self.test_lineplan.items()}
                for line in self.test_lineplan:
                    if line in self.linepool:
                        f[line].start = self.test_lineplan[line]
                        f[line].lb = self.test_lineplan[line]
                        f[line].ub = self.test_lineplan[line]

            # optional: test a schedule by fixing variables
            if self.test_schedule is not None:
                logging.info('Fixing test schedule variables')
                for event in self.test_schedule:
                    if event in Events:
                        pi[event].start = self.test_schedule[event]
                        pi[event].lb = self.test_schedule[event]
                        pi[event].ub = self.test_schedule[event]

            # optional: suggest gurobi solution for pi
            if self.start_pi is not None:
                for key, value in self.start_pi.items():
                    pi[key].start = value

            # optional: suggest gurobi solution for z
            if self.start_z is not None:
                for key, value in self.start_z.items():
                    z[key].start = value

            # optional: suggest gurobi solution for x
            if self.start_x is not None:
                for key, value in self.start_x.items():
                    x[key].start = value

            # optional: suggest gurobi solution for y
            if self.start_y is not None:
                for key, value in self.start_y.items():
                    y[key].start = value

            # objective function
            logging.info('Defining objective function')

            # COLPP part
            objective_line_planning = gp.quicksum(
                f[line] * self.linepool[line]['cost'] for line in self.linepool)

            # PESP part
            objective_time_tabling = 0
            for activity, properties in self.ean.Activities.items():
                weight = properties['weight']
                objective_time_tabling += weight * x[activity]

            # concatenate both parts
            model.setObjective(
                objective_line_planning + objective_time_tabling, gp.GRB.MINIMIZE)

            # Adding constraints
            logging.info('Adding constraints')

            # COLPP constraints
            edges_of_line = {edge: [] for edge in E}
            for line in self.linepool:
                stations = self.linepool[line]['stations']
                edges_between_stations = zip(
                    stations, stations[1:]
                )
                for edge in edges_between_stations:
                    edges_of_line[edge].append(f[line])

            # Now add for each edge the constraint the minimal and maximal frequency
            for edge in edges_of_line:
                model.addConstr(
                    gp.quicksum(edges_of_line[edge]) >= E[edge]['fmin'],
                    name=f'fmin_{edge}'
                )
                model.addConstr(
                    gp.quicksum(edges_of_line[edge]) <= E[edge]['fmax'],
                    name=f'fmax_{edge}'
                )
            del edges_of_line

            # time tabling constraints
            for activity, properties in self.ean.Activities.items():
                lb, ub = properties['lb'], properties['ub']
                i, j = activity[0], activity[1]
                model.addConstr(
                    pi[j] - pi[i] + self.T * z[activity] >= lb,
                    name=f'lb_{(i, j)}')
                model.addConstr(
                    pi[j] - pi[i] + self.T * z[activity] <= ub,
                    name=f'ub_{(i, j)}')

            # y - connection constraints
            for activity in self.ean.Activities: #y's
                line1, line2 = activity[0][2], activity[1][2]
                model.addConstr(
                    y[activity] <= f[line1],
                    name=f'y_{line1}_{activity}')
                model.addConstr(
                    y[activity] <= f[line2],
                    name=f'y_{line2}_{activity}')
                model.addConstr(
                    f[line1] + f[line2] <= y[activity] + 1,
                    name=f'y_{line1}_{line2}_{activity}')
                
            # x - connection constraint
            for activity, properties in self.ean.Activities.items():
                lb, ub = properties['lb'], properties['ub']
                i, j = activity[0], activity[1]
                model.addConstr(
                        x[activity] >= pi[j] - pi[i] + self.T*z[activity] + self.M*(y[activity] - 1),
                        name=f'x_{i}_{j}'
                    )
                
            # setting optimization parameters
            if self.TIMELIMIT is not None:
                model.setParam('TimeLimit', self.TIMELIMIT)

            # Initialize callback data storage
            callback_data = []
            start_time = time_module.time()
            model._callback = my_callback

            # optimize the model
            logging.info('Solving the model')
            model.optimize(my_callback)

            # read out optimal solution
            if model.status == gp.GRB.INFEASIBLE:
                logging.warning('Model is infeasible!')
                logging.info('Calculating smallest set of infeasible constraints')
                model.computeIIS()
                model.write(f'model_{model.modelName}.ilp')

            elif model.status == gp.GRB.OPTIMAL or model.status == gp.GRB.TIME_LIMIT:
                logging.info(f'Model terminated with the status {model.status}!')
                self.obj_val = model.objVal
                self.lineplan = {(line): (f[line].x) for line in self.linepool}
                self.schedule = {event: pi[event].x for event in Events}
                self.z_modulo_param = {
                    activity: z[activity].x for activity in Activities}
                self.y_connection_param = {
                    activity: y[activity].x for activity in Activities}
                self.x_connection_param = {
                    activity: x[activity].x for activity in Activities
                }

            # Save the callback data to a CSV file
            if self.SAVE_CALLBACK:
                with open('callback_data.csv', 'w', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(['Time', 'Objective Value'])
                    writer.writerows(callback_data)

            return
    

class LinTimApprox:
    '''
    A metamodel which solves the integrated line planning and timetabling
    problem.
    '''
    def __init__(
            self,
            ptn: PTN,
            ean: EAN,
            linepool: dict,
            T=60,
            PATIENCE=5,
            CSV_NAME='lintim_approx.csv',
            TIMELIMIT=None,
            MODE_TIME_INC=False,
            TIME_INCREASE=300,
            MAXTIME=3600) -> None:

        logging.info('INIT LinTim Approx')

        # input variables
        self.ptn = ptn
        self.ean = ean
        self.linepool = linepool
        self.T = T
        self.MAX_CONS_REJECTS = PATIENCE
        self.CSV_NAME = CSV_NAME

        # metamodel model variables
        self.beta = None
        self.beta_opt = None
        self.consecutive_rejects = 0
        self.i = 1
        self.start_time = time_module.time()
        self.simulation = LintimSolver(
            ptn = None,
            ean = None,
            linepool = None,
            T = self.T,
            run_solve = False
        )

        # output vars
        self.df = pd.DataFrame({})
        self.OPT_lineplan = None
        self.OPT_timetable = None
        self.objValMM = np.nan
        self.OPT_obj_val = np.infty

        # optimization params
        self.TIMELIMIT = TIMELIMIT
        self.MODE_TIME_INC = MODE_TIME_INC  # mode, that increases timelimit if we do not improve, bool 
        self.START_TIMELIMIT = TIMELIMIT  # reset timelimit to starting timelimit if we have improved in iteration
        self.TIME_INCREASE = TIME_INCREASE  # if we do not improve, increase timelimit by amount
        self.MAXTIME = MAXTIME  # upper bound for increasing timelimit

        # run the meta model
        self.run()

    def update_df(self, lintim_obj_val: float) -> None:
        '''
        Updates the stats of the metamodel which are stored in a
            dataframe.
        '''
        # creating new row
        data = {}
        data['nu(MM)'] = self.objValMM
        data['nu(LinTim)'] = lintim_obj_val
        data['nu*(LinTim)'] = self.OPT_obj_val
        data['sum(betas)'] = self.beta.sum()
        data['beta_pd'] = is_positive_definite(self.beta)
        data['time_finished'] = time_module.time() - self.start_time

        # updating dataframe and saving it
        self.df = pd.concat([self.df, pd.DataFrame([data])], ignore_index=True)
        self.df.to_csv(self.CSV_NAME)
        return

    def run(self) -> None:
        # INIT: betas
        logging.info('INIT: betas')

        # variant 1: take the identity matrix
        self.beta = np.eye(len(self.linepool))

        # variant 2: sum of weights of activities of intersecting lines
        # self.beta = pesp_weights_to_betas(self.ean, self.linepool, self.T)
        # self.beta = self.beta + np.eye(self.beta.shape[0], self.beta.shape[0])  # for p.d.

        while self.consecutive_rejects < self.MAX_CONS_REJECTS:
            logging.info(f'**Iteration {self.i}')

            # calculating COLPP
            logging.info('Calculating COLPP')
            stime = time_module.time()
            lpp_solver = ColppSolver_beta(
                ptn=self.ptn,
                linepool=self.linepool,
                beta=self.beta,
                run_solve=True,
                TIMELIMIT=self.TIMELIMIT,
            )
            self.objValMM = lpp_solver.obj_val
            logging.info(f'TIME_INFO MM: It took {time_module.time() - stime}')

            # simulating lineplan
            logging.info('Simulating lineplan')
            stime = time_module.time()
            lintim_solver = LintimSolver(
                ptn=self.ptn,
                ean=self.ean,
                linepool=self.linepool,
                T=self.T,
                test_lineplan=lpp_solver.lineplan,
                start_pi=self.simulation.schedule,
                start_z=self.simulation.z_modulo_param,
                start_x=self.simulation.x_connection_param,
                start_y=self.simulation.y_connection_param,
                TIMELIMIT=self.TIMELIMIT,
            )
            self.simulation = lintim_solver
            logging.info(f'TIME_INFO LinTim: It took {time_module.time() - stime} in iteration {self.i}')

            # calculating optimal betas
            logging.info('Calculating optimal betas')
            self.beta_opt = pesp_part_to_betas(
                ean=self.ean,
                schedule=lintim_solver.schedule,
                linepool=self.linepool,
                y_connection_param=lintim_solver.y_connection_param,
                z_modulo_param=lintim_solver.z_modulo_param
            )

            # updating betas
            logging.info('Updating betas')
            self.update_beta()

            logging.info('Checking if improved')
            if lintim_solver.obj_val < self.OPT_obj_val:
                # overwrite optimal stats and reset rejection counter
                logging.info('Objective value has been improved!')
                logging.info('Saving lp, tt and its objective value')
                self.OPT_lineplan = lintim_solver.lineplan
                self.OPT_timetable = lintim_solver.schedule
                self.OPT_obj_val = lintim_solver.obj_val
                self.consecutive_rejects = 0
                if self.MODE_TIME_INC:
                    self.TIMELIMIT = self.START_TIMELIMIT

            else:
                logging.info('Objective value has NOT been improved.')
                self.consecutive_rejects += 1
                if self.MODE_TIME_INC:
                    self.TIMELIMIT = min(
                        self.TIMELIMIT + self.TIME_INCREASE,
                        self.MAXTIME)

            # update and saving dataframe
            logging.info('Update and saving dataframe')
            self.update_df(lintim_solver.obj_val)

            # updating iteration counterS
            self.i += 1

        return
    
    def update_beta(self) -> None:
        # use entrywhise maximum
        self.beta = np.maximum(self.beta, self.beta_opt)

        return
