# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 02:50:12 2022  
@author: Rahman Khorramfar

This is the main script to run the optimization model.
"""

# from IPython import get_ipython;
# get_ipython().magic('reset -f') # to clear the namespace
# get_ipython().magic('clear');


import time;
import gurobipy as gp;
from gurobipy import GRB, quicksum,LinExpr;
from Setting import Setting;
from SystemClasses import EV,GV;
import numpy as np;
import sys; 

import multiprocessing as mp

#from Modules import Create_Power_System_Model;
#import Modules;
#%% Set Default Setting for the Porblem 
Setting.Power_network_size = 17; 
Setting.base_year = 2020; 
Setting.dispatch_year = 2020;
Setting.rep_day_folder = 'Representative_Days';#extreme_days
Setting.num_rep_days = 365;
Setting.emis_case = 4; #case 6, no biogas
Setting.electrification_scenario = 'Historical';   # HE, HX, ME, MX, RF, and Historical
Setting.emis_reduc_goal = 0;  
Setting.VRE_share = 0.0;
Setting.solver_gap = 0.01;
Setting.wall_clock_time_lim = 10; #hour
Setting.UC_active = 1;
Setting.relax_UC_vars = True;
Setting.relax_int_vars = 1;
Setting.solver_thread_num = 4;
Setting.e_shed_penalty = 1e4;
Setting.g_shed_penalty = 1e4;

Setting.Metal_air_storage_cost='no-metal-air';# low and high, no-metal-air
Setting.CCS_allowed = 1;
Setting.LCF_supply_curve = 0;
Setting.methane_leakage = 0;

if len(sys.argv)>1:
    print(str(sys.argv));
    Setting.rep_day_folder = sys.argv[1];
    Setting.Power_network_size = int(sys.argv[2]);
    Setting.base_year = int(sys.argv[3]);
    Setting.dispatch_year = int(sys.argv[4]);
    Setting.num_rep_days = int(sys.argv[5]);
    Setting.emis_case = int(sys.argv[6]);
    Setting.electrification_scenario = sys.argv[7];
    Setting.emis_reduc_goal = float(sys.argv[8]); # %80
    Setting.VRE_share = float(sys.argv[9]);
    Setting.solver_gap = float(sys.argv[10]);
    Setting.wall_clock_time_lim = int(sys.argv[11]); # hour
    Setting.UC_active = bool(int(sys.argv[12]));
    Setting.relax_UC_vars = bool(int(sys.argv[13]));
    Setting.relax_int_vars = bool(int(sys.argv[14]));
    Setting.solver_thread_num = int(sys.argv[15]);
    Setting.CCS_allowed = int(sys.argv[17]);
    
Setting.expansion_allowed = 0;
Setting.LCF_emissions = 0;
Setting.CRM_reserve = -1; # as recommended by NERC
Setting.wall_clock_time_lim = Setting.wall_clock_time_lim*3600; # convert to second for Gurobi
Setting.print_result_header = 0;
Setting.copper_plate_approx = 0; 
Setting.print_all_vars = 0;
s_time = time.time();

#%% recall the heuristic to obtain an approximate solution

import Module;



    
def run_optimization(scenario_name):
    # try:
    # Set the scenario name
    Setting.disruption_scenario = scenario_name
    # if scenario_name == 'Scen_0':
    #     Setting.print_result_header = 1
    # else:
    #     Setting.print_result_header = 0
    # 
    #%% create and recall and run the modules
    Model = gp.Model()
    Module.Power_System_Model(Model); # Power System
    Module.NG_System_Model(Model);# NG System
    Module.Coupling_constraints(Model);# Coupling Constraints


    #% add objective function and run
    Model.modelSense = GRB.MINIMIZE;
    #Model.setObjective(EV.e_system_cost);
    if Setting.emis_case==1:
        Model.setObjective(EV.e_system_cost);
    else:
        Model.setObjective(GV.g_system_cost+ EV.e_system_cost);
    Model.setParam('OutputFlag', 0);
    # Model.setParam('MIPGap', Setting.solver_gap);
    Model.setParam('Timelimit', Setting.wall_clock_time_lim);

    # Model.setParam('Threads',Setting.solver_thread_num);
    #Model.setParam('Presolve',2); # -1 to 2
    # Model.setParam('MIPFocus',2); # 0 to 3
    # Model.setParam('NodefileStart',0.5);
    # Model.setParam('PreSparsify',1); 

    # After setting up the model but before Model.optimize()
    Model.update()  # Finalizes any pending model modifications


    t5 = time.time()
    Model.optimize();
    t6 = time.time()
    print(f'objective value: {np.round(Model.ObjVal/1e9,2)}e9');
    print(f"Num of Vars: {Model.NumVars}");
    print(f"Num of Int Vars: {Model.NumIntVars}");
    print(f"Num of Constraints: {Model.NumConstrs}");

    MIP_gap=0;
    # if Setting.relax_int_vars==False:
    # MIP_gap = 100*Model.MIPGap;

    #% print
    Module.Get_var_vals(Model);
    Module.Publish_results(s_time,MIP_gap);

    print(sys.argv);
    print(f"\n\n Elapsed time (seconds): {time.time()-s_time}");
    print(f"MIP Gap (%): {MIP_gap}");
    # print(f"\n Establishment cost: {format(EV.est_cost_val,'.2E')}");
    # print(f"Decommissioning cost: {format(EV.decom_cost_val,'.2E')}");
    # print(f"FOM cost: {format(EV.FOM_cost_val,'.2E')}");
    # print(f"VOM cost: {format(EV.VOM_cost_val,'.2E')}");
    # print(f"Fuel cost: {format(EV.nuc_fuel_cost_val,'.2E')}");
    # print(f"Fuel cost: {format(EV.gas_fuel_cost_val,'.2E')}");
    # print(f"Startup cost: {format(EV.startup_cost_val,'.2E')}");
    # print(f"Sheding cost: {format(EV.shedding_cost_val,'.2E')}");
    # print(f"Storage cost: {format(EV.elec_storage_cost_val,'.2E')}");
    # print(f"Trans. Est. cost: {format(EV.est_trans_cost_val,'.2E')}");
    # print(f"e emission: {format(EV.emis_amount_val,'.2E')}");
    # print(f"CCS Cost: {format(EV.CCS_cost_val,'.2E')}");
    # print(f" Power system cost: {format(EV.e_system_cost_val,'.3E')}\n")

    # print(f"\n\n storage inv. cost: {format(GV.inv_str_cost_val,'.2E')}");
    # print(f"pipeline inv. cost: {format(GV.inv_pipe_cost_val,'.2E')}");
    # print(f"Shedding cost: {format(GV.shed_cost_val,'.2E')}")
    # print(f"RNG cost: {format(GV.RNG_cost_val,'.2E')}")
    # print(f"FOM storage cost: {format(GV.fom_str_cost_val,'.2E')}")
    # print(f"NG import cost: {format(GV.import_cost_val,'.2E')}")
    # print(f"NG Emission: {format(GV.emis_amount_val,'.2E')}")
    # print(f"NG system cost: {format(GV.g_system_cost_val,'.3E')}")
    print(f"Total Energy cost: {format(Model.ObjVal,'.3E')}")
    print(f"best lower bound: {format(Model.ObjBound,'.3E')}")

    # except Exception as e:
    #     print(f"Error in scenario {scenario_name}: {e}")


def parallel_optimization(scenario_list):
    start_time = time.time()

    # Create a pool of processes, with one process per CPU core (or modify as needed)
    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.map(run_optimization, scenario_list)

    print(f"Total execution time: {time.time() - start_time:.2f} seconds")



if __name__ == "__main__":
    # Define the list of scenarios you want to run
    #scenarios = ['Scen_0']
    scenarios = []
    for n in range(1, 10):
        for i in range(1, 2):
            scenarios.append(f'Trans_scenario__n_{2*n}_#{i}')

    scenarios = ['Trans_scenario__n_1_#1', 'Trans_scenario__n_10_#1']

    run_optimization('Scen_0')
    # Run parallel optimization
    #parallel_optimization(scenarios)