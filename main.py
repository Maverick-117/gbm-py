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
compDosesQ = True; #true: yes you're comparing doses; false: no, you're comparing reprogramming versus non-reprogramming
kimICQ = True; # true: use kim's IC; false: use Yu's IC
c_dep_sQ = False; # true: make c in radiotherapy dependent on s; false: don't do that
kimReprogQ = False; 
kimDeathValQ = True; # true: sets death value to what Kim used; false: sets death value equal to 
deathFdbkQ = True; # false: set death feedback gain to 0; true: don't use the false option
    
model0Q = True; # true: using vo's no-survivin model and equations; false: not using it 
model1Q = False; # true: using vo's survivin-protects model and equations; false: not using it
model2Q = False; # true: using vo's survivin-protects-dedifferentiates model and equations; false: not using it

### VARS RELATED TO EXPERIMENTAL QUESTS ###
goldenLaptopQ = False;


# Radiotherapy model
a,b =  np.array([0.17, 0.02]); 
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
h1 = 10**5; # feedback on css div 
h2 = 10**5; # feedback on dcc div
pwr = 3;#Inf;
ROI_radius = 1; # Radius of the simulation region of intrest
rho = 10**9; # density of cells in the region of interest (cells/cm^3)
total_cell_num = 4/3*np.pi*(ROI_radius ** 3)*rho;
post_therapy_end = 1100;
#srv_start = 0;
mu_start = 0;
ss = [6];#np.arange(12).tolist();#[6];#[4,8,9,5,10,11];
#[0,1,3,4,8,9,5,10,11];#[4,8,9,5,10,11];[0,1,2,3,6,7]
z = 1; n = 1;
l_vec =  [0,  0, l_w, l_s, l_w , l_s , 0 ,0 ,l_w,l_w,l_s,l_s];
h1_vec = [0, h1, 0  , 0  , h1  , h1  , h1,0 ,h1 ,0  , h1, 0 ];
h2_vec = [0, h2, 0  , 0  , h2  , h2  , 0 ,h2,0  ,h2 , 0 , h2];

# man-behind-the-curtain setup.
switch_vec = [subSelectQ, srvQ, compDosesQ, deathFdbkQ, c_dep_sQ, kimReprogQ, kimDeathValQ, kimICQ, model0Q, model1Q, model2Q];
misc_pars = [DT, post_therapy_end, pwr, h1, h2, l_vec, ss, a, b];
par_setup_vec, string_setup_vec = funciones.parameter_setup(switch_vec, misc_pars);
#total_start_frac, d, Doses, Frac, C, cont_p_a, chi, mu_bar, hd, zeta_mult1, zeta_mult2, cont_c, rng = par_setup_vec;
[total_start_frac, d, Doses, Frac, C, cont_p_a, chi, mu_bar, hd, cont_c, rng, mu_start] = par_setup_vec;
hd_str_mod, v_suffix, v_suffix_save_A, v_suffix_save, fg_subtitle, reprog_suffix, hd_str, hd_suffix, group_name, c_fdbk_mod = string_setup_vec;

#stray parameters to setup
total_start_cell_num = total_cell_num*total_start_frac;
cont_p_b = 10*cont_p_a; compt_mult = 100; # control the strength of inhibitory signal
#srvn_zeta = [3.6 * zeta_mult1, 0.05 * zeta_mult2]; # original values are empirically derived
srvn_zeta = [3.6, 0.05]; # original values are empirically derived
srvn_csc = srvn_zeta[0]; srvn_dcc = srvn_zeta[1];
surv_vec = [cont_p_a, cont_p_b, compt_mult, srvn_csc, srvn_dcc, cont_c]; #assuming these control parameters are constant
# surv_vec_list = [];
# zetas = [[3.6,0.05]];#[[0,0],[3.6,0.05],[0.05, 3.6]]
# for ar in range(len(zetas)):
#     surv_vec_list+=surv_vec[:3]+zetas[ar]+[surv_vec[-1]];
    
# stray strings to setup
fdbk_type_vec = ['No', 'Div Only', 'Weak', 'Strong', 'Weak', 'Strong'];
title_vec = [" nothing", " m_u, m_v", " p (weak)", " p (strong)", " m_u, m_v" + " and " + " p (weak)", " m_u, m_v" + " and " + " p (strong)"] + [hd_suffix];
color = ['k','r','b','m'];
cell_lines = ["U373MG"];

#### VARIATION ZONE FOR PARAMETERS
d *= 1;
time_pts1 = 200;
time_pts2 = 200;
day_month = "31_Oct"; #05_July last used
base_model_name = 'k2_model'
model_suffix = "conventional_BED\\"#"comparing_conventionals\\";
case = "reprod_old2";#"All\\";
if goldenLaptopQ:
    base_dirGD = "C:\\Users\\jhvo9\\Google Drive (vojh1@uci.edu)"
else:
    base_dirGD = "G:\\My Drive"
drty = base_dirGD + "\\a PhD Projects\\GBM Modeling\\python scripts\\data\\"+base_model_name+"\\"+model_suffix+"\\"+case+"\\"+day_month; # _div_rate_diff

if kimDeathValQ:
    deathVal_dir = '\\death_val_of_kim';
else:
    deathVal_dir = '\\death_val_of_not_kim';
if kimReprogQ:
    sub_drty = "\\kim_reprog";
else:
    sub_drty = "\\corrected_reprog";
total_drty = drty + deathVal_dir + sub_drty;
if not exists(total_drty):
    makedirs(total_drty);
if compDosesQ:
    comp_str = "Gy";
    comp_list = Doses;
    Frac_list = Frac;
else:
    comp_str = "reprog";
    comp_list = C;
    Frac_list = [Frac, Frac];
###
if deathFdbkQ:
    deathFdbk_str="_w_death_fdbk";
else:
    deathFdbk_str="_w_no_death_fdbk";
print(deathVal_dir + '\\' + sub_drty+comp_str+deathFdbk_str);
if compDosesQ:
    variegator = Doses; #outer loop will put the unchanging vector first
    fixor = C; # inner loop will use what's changing
    # technically this is unnecessary for data generation
    # but for plotting this helps
    # this could be easily generalized for any C and Doses used
else:
    variegator = C;
    fixor = Doses;
case_length = len(variegator);

for lll in rng:
    print("index val:", lll,"\n")
    u_sc = list(range(case_length)); u_dc = list(range(case_length)); u_srv = list(range(case_length)); t_vec = list(range(case_length));
    un_sc = list(range(case_length)); un_dc = list(range(case_length)); un_srv = list(range(case_length)); tn_vec = list(range(case_length));
    l = l_vec[lll]; h1 = h1_vec[lll]; h2 = h2_vec[lll];
    ## Tumor Growth ODE and radiotherapy simulation
    # for vv, el in surv_vec_list:
        # surv_vec_list = el;
    for gg in range(len(cell_lines)):
        print('Cell Line: '+cell_lines[gg]+'\n')
        print('a = '+str(a)+' b = '+str(b)+'\n');
        for g in range(len(fixor)):
            if compDosesQ:
                fix_str = str(C[g])+"_reprog";
            else:
                fix_str = str(Doses[g])+"_Gy";
                
            print("g = ",g)
            if compDosesQ:
                c = C[g];
            else:
                D = Doses[g];
                frac_num = Frac[g];
            for k in range(case_length):
                print("k = ",k)
                ## looping over reprogramming values
                 # Defining Variables
                if compDosesQ:
                    D = Doses[k];
                    frac_num = int(Frac[k]);
                else:
                    c = C[k];
                    
                sc_start = total_start_frac*F;
                tc_start = total_start_frac-sc_start;
                total_cell_num = 4/3*np.pi*10**9
                total_start_num = total_start_frac*total_cell_num
                
                # Defining treatment days and ODE simulation start time after each fraction
                weeks = int(frac_num//5);
                total_days = frac_num + 2*weeks;
                acq_end = treat_start + post_therapy_end;#treat_start + total_days + acq_days_after_RT - 1;
                treat_days = np.tile([1,1,1,1,1,0,0], weeks+1)[:total_days].nonzero()[0]+treat_start;
                sim_resume_days = treat_days+10/(60*24); # ODE simulation resumes 10 minutes after RT
                treat_days = np.array(treat_days.tolist() + [acq_end]);  
                LQ_param = [a1, b1, a2, b2, c, D];
                
                sim_values = [model0Q, model1Q, model2Q, kimReprogQ, total_cell_num, treat_days, mu_start, LQ_param, total_start_frac, sc_start, sim_resume_days, surv_vec, time_pts1, time_pts2];
                para_values = [r1, r2, d, p, h1, h2, hd, z, l, n, sig, mu_bar, chi];
                
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
                    print("Saving", total_drty+"\\TU_"+fix_str+str(comp_list[k])+'_'+comp_str+'_'+str(Frac_list[k])+'_Fracs_'+str(post_therapy_end)+'days_'+deathFdbk_str+"_"+str(lll)+"nn.txt")
                    np.savetxt(total_drty+"\\TU_"+fix_str+str(comp_list[k])+'_'+comp_str+'_'+str(Frac_list[k])+'_Fracs_'+str(post_therapy_end)+'days_'+deathFdbk_str+"_"+str(lll)+"nn.txt",[t_vec[k], u_sc[k], u_dc[k], u_srv[k]],header="Time, CSC, DCC, SRV")
                    print("Saving", total_drty+"\\TU_none_"+fix_str+str(comp_list[k])+'_'+comp_str+'_'+str(Frac_list[k])+'_Fracs_'+str(post_therapy_end)+'days_'+deathFdbk_str+"_"+str(lll)+"nn.txt")
                    np.savetxt(total_drty+"\\TU_none_"+fix_str+str(comp_list[k])+'_'+comp_str+'_'+str(Frac_list[k])+'_Fracs_'+str(post_therapy_end)+'days_'+deathFdbk_str+"_"+str(lll)+"nn.txt",[tn_vec[k], un_sc[k], un_dc[k], un_srv[k]],header="Time_n, CSC_n, DCC_n, SRV_n")
