import os
import pandas as pd
import numpy as np
from Setting import Setting
#from ProblemData import time_weights

DATA_path = os.path.dirname(os.getcwd())+"/DATA"

def time_weights(nD, base_year):
    """Loads the representative day data
        INPUTS :
            nD (int): number of representative days
            base_year (int): year of the modeling
        OUTPUTS :
            e_time_weight (array): time weights for each hour of the representative days (used for electricity)
            g_time_weight (array):  time weights for each representative day (used for gas)
            g_rep_days (array):   indices of the representative days (used for gas)
            e_rep_hrs (array): representative hours (used for electricity)
            days_in_cluster (dict):  mapping cluster index to days in that cluster
            days2Medoid (dict):  mapping each day to its corresponding medoid
    """ 
    #df = pd.read_csv(os.getcwd()+r'\Rep_Days\num_rep_days='+str(nD)+'.csv');
    #df = pd.read_csv(os.getcwd() + '/' + fName+'/'+fName +'_medians_k='+str(nD)+'.csv');
    
    #df = pd.read_csv(os.getcwd() + '/' + fName+'/Rep-Days-Year='+str(Setting.base_year)+'/cluster_year='+str(Setting.base_year)+'_days='+str(nD)+'.csv');

    
    df = pd.read_csv(DATA_path + f'/Representative_Days/cluster_year={base_year}_days={nD}.csv')

    s1 = df[df['is_medoid?']==1]
    s2 = s1.sort_values(['Cluster'])

    g_rep_days = np.array(s2.index)
    e_time_weight = np.zeros(24*nD,dtype=int)
    e_rep_hrs = np.zeros(24*nD,dtype=int)
    g_time_weight = np.zeros(nD,dtype=int)

    days_in_cluster = dict()
    days2Medoid = dict()
    
    # sort g_rep_days and g_weights
    g_rep_days = sorted(g_rep_days)
    for i in range(nD):
        s1 = np.array(df['Cluster'])
        i2 = df['Cluster'].iloc[g_rep_days[i]]
        s2 = np.where(s1==i2)
        g_time_weight[i]= len(s2[0])
        for d in s2[0]:
            days2Medoid[d] = g_rep_days[i]
        days_in_cluster[i] = s2[0]
        e_time_weight[i*24:(i+1)*24] += g_time_weight[i]
        for j in range(24):
             e_rep_hrs[i*24+j] =g_rep_days[i]*24+j

    return e_time_weight, g_time_weight, g_rep_days,e_rep_hrs, days_in_cluster, days2Medoid

e_time_weight, g_time_weight, g_rep_days,e_rep_hrs, days_in_cluster, days2Medoid = time_weights(365, 2020)

print(days2Medoid)

