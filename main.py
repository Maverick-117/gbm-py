# -*- coding: utf-8 -*-
"""
Created on Fri May  7 19:52:20 2021

@author: jhvo9
"""

import numpy as np
import funciones
from os import makedirs
from os.path import exists
#### SETUP

saveQ = True;
srvQ = False; # true means that survivin is at play for any sim; false means that survivin is barred from any sim
subSelectQ = True; # true: downselect on a particular kind of feedback regime; false: use it all
compDosesQ = False; #true: yes you're comparing doses; false: no, you're comparing reprogramming versus non-reprogramming
kimICQ = True; # true: use kim's IC; false: use Yu's IC
c_dep_sQ = False; # true: make c in radiotherapy dependent on s; false: don't do that
kimReprogQ = False; 
kimDeathValQ = False; # true: sets death value to what Kim used; false: sets death value equal to 

model0Q = True; # true: using vo's no-survivin model and equations; false: not using it 
model1Q = False; # true: using vo's survivin-protects model and equations; false: not using it
model2Q = False; # true: using vo's survivin-protects-dedifferentiates model and equations; false: not using it

### VARS RELATED TO EXPERIMENTAL QUESTS ###

probFdbkQ  = False;
divR1FdbkQ = False;
divR2FdbkQ = False;
deathFdbkQ = False; # false: set death feedback gain to 0; true: don't use the false option


# Radiotherapy model
a1 = 0.01
b1 = 1.77e-7
a2 = 0.125
b2 = 0.028
F = 0.016
sig =  0.4750;
DT =  3.9; # doubling time in days in order labeled by varible "cell_lines"
treat_start = 100; # treatment start time relative to start of simulation
acq_days_after_RT = 110; # simulation length after last fraction

# DE parameters
r1 = np.log(2)/DT; # growth rate of CSC
r2 = np.log(2)/DT; # growth rate of DCC
p = .505; # self renewal probability 
l_w = 10**(-7); # weak feedback on prob
l_s = 10**3; # strong feedback on prob
h = 10**5; # feedback on div
pwr = 3;#Inf;
ROI_radius = 1; # Radius of the simulation region of intrest
rho = 10**9; # density of cells in the region of interest (cells/cm^3)
total_cell_num = 4/3*np.pi*(ROI_radius ** 3)*rho;
post_therapy_end = 1000;
srv_start = 0;
ss = 4;
z = 1; n = 1;
l_vec         = [0  , 0, l_w, l_s, l_w, l_s];
h_vec         = [0  , h, 0  , 0  , h  , h  ];

# man-behind-the-curtain setup.
switch_vec = [subSelectQ, srvQ, compDosesQ, deathFdbkQ, c_dep_sQ, kimReprogQ, kimDeathValQ, kimICQ, model0Q, model1Q, model2Q];
misc_pars = [DT, post_therapy_end, pwr, h, l_vec, ss];
par_setup_vec, string_setup_vec = funciones.parameter_setup(switch_vec, misc_pars);
total_start_frac, d, Doses, Frac, C, cont_p_a, chi, mu_bar, hd, zeta_mult1, zeta_mult2, cont_c, rng = par_setup_vec;
hd_str_mod, v_suffix, v_suffix_save_A, v_suffix_save, fg_subtitle, reprog_suffix, hd_str, hd_suffix, group_name, c_fdbk_mod = string_setup_vec;

#stray parameters to setup
total_start_cell_num = total_cell_num*total_start_frac;
cont_p_b = 10*cont_p_a; compt_mult = 100; # control the strength of inhibitory signal
srvn_zeta = [3.6 * zeta_mult1, 0.05 * zeta_mult2];
srvn_csc = srvn_zeta[0]; srvn_dcc = srvn_zeta[1];
surv_vec = [cont_p_a, cont_p_b, compt_mult, srvn_csc, srvn_dcc, cont_c]; #assuming these control parameters are constant

# stray strings to setup
fdbk_type_vec = ['No', 'Div Only', 'Weak', 'Strong', 'Weak', 'Strong'];
title_vec = [" nothing", " m_u, m_v", " p (weak)", " p (strong)", " m_u, m_v" + " and " + " p (weak)", " m_u, m_v" + " and " + " p (strong)"] + [hd_suffix];
color = ['k','r','b','m'];
cell_lines = ["U373MG"];

#### VARIATION ZONE FOR PARAMETERS
d *= 1;
time_pts1 = 200;
time_pts2 = 200;
day_month = "17_June";
model_suffix = "_positive_death_feedback_expts";
drty = "C:\\Users\\jhvo9\\Google Drive (vojh1@uci.edu)\\a PhD Projects\\GBM Modeling\\python scripts\\data\\kim_model"+model_suffix+"\\"+day_month; # _div_rate_diff


if probFdbkQ:
    probFdbk_dir = '\\with_prob_fdbk';
else:
    probFdbk_dir = '\\without_prob_fdbk';

if divR1FdbkQ:
    divR1Fdbk_dir = '\\with_divR1_fdbk';
else:
    divR1Fdbk_dir = '\\without_divR1_fdbk';
    

if divR2FdbkQ:
    divR2Fdbk_dir = '\\with_divR2_fdbk';
else:
    divR2Fdbk_dir = '\\without_divR2_fdbk';
fdbk_dirs = probFdbkQ + divR1FdbkQ + divR2FdbkQ;
if kimDeathValQ:
    deathVal_dir = '\\death_val_of_kim';
else:
    deathVal_dir = '\\death_val_of_not_kim';
if kimReprogQ:
    sub_drty = "\\kim_reprog";
else:
    sub_drty = "\\corrected_reprog";
total_drty = drty + deathVal_dir + sub_drty + fdbk_dirs;
if not exists(total_drty):
    makedirs(total_drty);
if compDosesQ:
    ending1 = str(Doses[0])+"Gy";
    ending2 = str(Doses[1])+"Gy";
else:
    ending1 = str(C[0])+"reprog";
    ending2 = str(C[1])+"reprog";
###
if deathFdbkQ:
    ending1+="_w_death_fdbk";
    ending2+="_w_death_fdbk";
else:
    ending1+="_w_no_death_fdbk";
    ending2+="_w_no_death_fdbk";
print(deathVal_dir + '\\' + sub_drty+ending2);
print(deathVal_dir + '\\' + sub_drty+ending1);
if compDosesQ:
    variegator = Doses; #outer loop will put the unchanging vector first
    fixor = C; # inner loop will use what's changing
    # technically this is unnecessary for data generation
    # but for plotting this helps
    # this could be easily generalized for any C and Doses used
else:
    variegator = C;
    fixor = Doses;

for lll in rng:
    print("index val:", lll,"\n")
    u_sc = [0,1]; u_dc = [0,1]; u_srv = [0,1]; t_vec = [0,1];
    un_sc = [0,1]; un_dc = [0,1]; un_srv = [0,1]; tn_vec = [0,1];
    l = l_vec[lll]; h = h_vec[lll];
    ## Tumor Growth ODE and radiotherapy simulation
    for gg in range(len(cell_lines)):
        a,b =  np.array([0.17, 0.02]); 
        print('Cell Line: '+cell_lines[gg]+'\n')
        print('a = '+str(a)+' b = '+str(b)+'\n');
        for g in range(len(fixor)):
            print("g = ",g)
            if compDosesQ:
                c = C[g];
            else:
                D = Doses[g];
                frac_num = Frac[g];
            for k in range(len(variegator)):
                print("k = ",k)
                ## looping over reprogramming values
                 # Defining Variables
                if compDosesQ:
                    D = Doses[k];
                    frac_num = Frac[k];
                else:
                    c = C[k];
                    
                sc_start = total_start_frac*F;
                tc_start = total_start_frac-sc_start;
                total_cell_num = 4/3*np.pi*10**9
                total_start_num = total_start_frac*total_cell_num
                
                # Defining treatment days and ODE simulation start time after each fraction
                weeks = frac_num//5;
                total_days = frac_num + 2*weeks;
                acq_end = treat_start + post_therapy_end;#treat_start + total_days + acq_days_after_RT - 1;
                treat_days = np.tile([1,1,1,1,1,0,0], weeks+1)[:total_days].nonzero()[0]+treat_start;
                sim_resume_days = treat_days+10/(60*24); # ODE simulation resumes 10 minutes after RT
                treat_days = np.array(treat_days.tolist() + [acq_end]);  
                LQ_param = [a1, b1, a2, b2, c, D];
                
                sim_values = [model0Q, model1Q, model2Q, kimReprogQ, total_cell_num, treat_days, srv_start, LQ_param, total_start_frac, sc_start, sim_resume_days, surv_vec, time_pts1, time_pts2];
                para_values = [r1, r2, d, p, 1e0 * h, 1e0 * h, hd, z, l, n, sig, mu_bar, chi];
                
                U, T, U_none, T_none = funciones.dynamics(para_values,sim_values);
                u_sc[k] = U[0,:];
                u_dc[k] = U[1,:];
                u_srv[k] = U[2,:];
                t_vec[k] = T;
                un_sc[k] = U_none[0,:];
                un_dc[k] = U_none[1,:];
                un_srv[k] = U_none[2,:];
                tn_vec[k] = T_none;

    #### Saving Data
    if saveQ:
        np.savetxt(total_drty+"\\TU_"+ending1+"_"+str(lll)+".txt",[t_vec[0], u_sc[0], u_dc[0], u_srv[0]],header="Time, CSC, DCC, SRV")
        np.savetxt(total_drty+"\\TU_"+ending2+"_"+str(lll)+".txt",[t_vec[1], u_sc[1], u_dc[1], u_srv[1]],header="Time, CSC, DCC, SRV")
        np.savetxt(total_drty+"\\TU_none_"+ending1+"_"+str(lll)+".txt",[tn_vec[0], un_sc[0], un_dc[0], un_srv[0]],header="Time_n, CSC_n, DCC_n, SRV_n")
        np.savetxt(total_drty+"\\TU_none_"+ending2+"_"+str(lll)+".txt",[tn_vec[1], un_sc[1], un_dc[1], un_srv[1]],header="Time_n, CSC_n, DCC_n, SRV_n")
