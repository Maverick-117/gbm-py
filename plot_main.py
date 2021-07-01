# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 11:11:18 2021

@author: jhvo9
"""
### TODO:
    # make axis font bigger
    # set some kind of axis limits?

import matplotlib.pyplot as plt;
import numpy as np;
from os import makedirs
from os.path import exists
total_cell_num = 4/3*np.pi*(1 ** 3)*1e9;
goldenLaptopQ = True;
compDosesQ = False; #true: yes you're comparing doses; false: no, you're comparing reprogramming versus non-reprogramming
saveQ = True;
log_yscaleQ = True;
r2 = np.log(2)/3.9;
Doses = [2,1]; n = 1;
C = [0, 0.005196];
subSelectQ = True; kimReprogQ = False; kimDeathValQ = True; 
l_w = 10**(-7); # weak feedback on prob
l_s = 10**1; # strong feedback on prob
h1 = 10**5; # feedback on css div 
h2 = 10**5; # feedback on dcc div
l_vec =  [0,  0, l_w, l_s, l_w , l_s , 0 ,0 ,l_w,l_w,l_s,l_s];
h1_vec = [0, h1, 0  , 0  , h1  , h1  , h1,0 ,h1 ,0  , h1, 0 ];
h2_vec = [0, h2, 0  , 0  , h2  , h2  , 0 ,h2,0  ,h2 , 0 , h2];
#pre_lll_load = 9;
if subSelectQ:
    #lll_vec = [pre_lll_load];
    #lll_vec = [4,8,9];
    lll_vec = [0,1,3,4,8,9,5,10,11];#[4,8,9,5,10,11];
else:
    lll_vec = list(range(len(l_vec)));
if log_yscaleQ:
    log_str = " (log)";
    log_dir = "\\log_yscale\\";
else:
    log_str = "";
    log_dir = "\\linear_yscale\\"
### TO VARY
deathFdbkQ = True; # false: set death feedback gain to 0; true: don't use the false option

### SETTINGS
date_dir = '23_June\\'; 
model_suffix = "_death_modification_2";

if kimDeathValQ:
    deathVal_dir = 'death_val_of_kim\\';
else:
    deathVal_dir = 'death_val_of_not_kim\\';

if kimReprogQ:
    sub_drty = "kim_reprog\\";
else:
    sub_drty = "corrected_reprog\\";
if compDosesQ:
    unit = ' Gy';
    csr = Doses;
else:
    unit = ' reprog'
    csr = C;
color = ['k','r','g','b','m'];
if compDosesQ:
    ending1 = str(Doses[0])+unit[1:];
    ending2 = str(Doses[1])+unit[1:];
else:
    ending2 = str(C[0])+unit[1:];
    ending1 = str(C[1])+unit[1:];

    
if deathFdbkQ:
    hd = 32520.32520325203;#1e5; 
    t1 = "_w_death_fdbk";
    ending1+=t1;
    ending2+=t1;
else:
    hd = 0.0; 
    t1= "_w_no_death_fdbk";
    ending1+=t1;
    ending2+=t1;

### PLOTTING
for lll_load in lll_vec:
    print("evaluating and graphic index:", lll_load)
    l = float(l_vec[lll_load]);
    h1 = float(h1_vec[lll_load]);
    h2 = float(h2_vec[lll_load]);
    plot_str2, plot_str1 = "TU_"+ending2+"_"+str(lll_load), "TU_"+ending1+"_"+str(lll_load);
    plotn_str2, plotn_str1 = "TU_none_"+ending2+"_"+str(lll_load), "TU_none_"+ending1+"_"+str(lll_load);
    probFdbkQ = l > 0;
    divR1FdbkQ = h1 > 0;
    divR2FdbkQ = h2 > 0;
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
    fdbk_dirs = probFdbk_dir + divR1Fdbk_dir + divR2Fdbk_dir;
    if goldenLaptopQ:
        base_dirGD = "C:\\Users\\jhvo9\\Google Drive (vojh1@uci.edu)"
    else:
        base_dirGD = "G:\\My Drive"
    data_drty = base_dirGD+"\\a PhD Projects\\GBM Modeling\\python scripts\\data\\kim_model"+model_suffix+"\\"+date_dir+ deathVal_dir + sub_drty
    plot_drty = base_dirGD+"\\a PhD Projects\\GBM Modeling\\python scripts\\plots\\kim_model"+model_suffix+"\\"+date_dir+deathVal_dir + sub_drty + fdbk_dirs + log_dir
    if not exists(plot_drty):
        makedirs(plot_drty);
    TU_2 = np.loadtxt(data_drty+plot_str2+".txt");
    TU_1 = np.loadtxt(data_drty+plot_str1+".txt");
    TU_2_none = np.loadtxt(data_drty+plotn_str2+".txt");
    TU_1_none = np.loadtxt(data_drty+plotn_str1+".txt");
    
    t_vec = [TU_2[0,:], TU_1[0,:]];
    u_sc = [TU_2[1,:], TU_1[1,:]];
    u_dc = [TU_2[2,:], TU_1[2,:]];
    u_srv = [TU_2[3,:], TU_1[3,:]];
    
    tn_vec = [TU_2_none[0,:], TU_1_none[0,:]];
    un_sc = [TU_2_none[1,:], TU_1_none[1,:]];
    un_dc = [TU_2_none[2,:], TU_1_none[2,:]];
    un_srv = [TU_2_none[3,:], TU_1_none[3,:]];
    
    #T_stop = 200;
    T_stop = int(t_vec[0].max());
    
    st2_stop = np.flatnonzero(t_vec[0] < T_stop).max();
    st1_stop = np.flatnonzero(t_vec[1] < T_stop).max();
    st0_stop = np.flatnonzero(tn_vec[0] < T_stop).max();
    
    s_title = "(l:"+f"{l:.3}"+";h1:"+f"{h1:.3}"+";h2:"+f"{h2:.3}"+";hd:"+f"{hd:.3}"+")"
    plot_save_suffix = '_'+str(T_stop)+'days_'+t1+'_'+str(lll_load);
    
    figT, axT = plt.subplots(figsize=(8,8));
    if log_yscaleQ:
        axT.set_yscale('log')
    axT.plot(t_vec[0][:st2_stop],u_sc[0][:st2_stop] * total_cell_num,color[0]+'--',label='CSC '+str(csr[0])+unit)
    axT.plot(t_vec[0][:st2_stop],u_dc[0][:st2_stop] * total_cell_num,color[0]+'-',label='DCC '+str(csr[0])+unit)
    axT.plot(t_vec[1][:st1_stop],u_sc[1][:st1_stop] * total_cell_num,color[1]+'--',label='CSC '+str(csr[1])+unit)
    axT.plot(t_vec[1][:st1_stop],u_dc[1][:st1_stop] * total_cell_num,color[1]+'-',label='DCC '+str(csr[1])+unit)
    axT.plot(tn_vec[0][:st0_stop],un_sc[0][:st0_stop] * total_cell_num,color[2]+'--',label='CSC, no treatment')
    axT.plot(tn_vec[0][:st0_stop],un_dc[0][:st0_stop] * total_cell_num,color[2]+'-',label='DCC, no treatment')
    axT.legend();
    axT.set_ylabel("Cell Number"+log_str)
    axT.set_title(s_title)
    if saveQ:
        figT.savefig(plot_drty+"\\pop_plot_"+plot_save_suffix+".png",dpi=300)
    
    # TODO: maybe add the other ones too?
    csc_p = [];
    csc_p += [u_sc[0]/(u_sc[0] + u_dc[0])];
    csc_p += [u_sc[1]/(u_sc[1] + u_dc[1])];
    csc_p += [un_sc[0]/(un_sc[0] + un_dc[0])];
    figC, axC = plt.subplots(figsize=(8,8));
    if log_yscaleQ:
        axC.set_yscale('log')
    axC.plot(t_vec[0][:st2_stop],csc_p[0][:st2_stop],color[0]+'-',label='CSC '+str(csr[0])+unit)
    axC.plot(t_vec[1][:st1_stop],csc_p[1][:st1_stop],color[1]+'-',label='CSC '+str(csr[1])+unit)
    axC.plot(tn_vec[0][:st0_stop],csc_p[2][:st0_stop],color[2]+'-',label='CSC, no treatment')
    axC.legend();
    axC.set_ylabel("CSC frac")
    axC.set_title(s_title)
    if saveQ:
        figC.savefig(plot_drty+"\\csc_plot_"+plot_save_suffix+".png",dpi=300)
    
    figD, axD = plt.subplots(figsize=(8,8));
    if log_yscaleQ:
        axD.set_yscale('log')
    axD.plot(t_vec[0][:st2_stop],hd*u_dc[0][:st2_stop]**n/(1+hd*u_dc[0][:st2_stop]**n),color[0]+'-',label='death feedback '+str(csr[0])+unit)
    axD.plot(t_vec[1][:st1_stop],hd*u_dc[1][:st1_stop]**n/(1+hd*u_dc[1][:st1_stop]**n),color[1]+'-',label='death feedback '+str(csr[1])+unit)
    axD.plot(tn_vec[0][:st0_stop],hd*un_dc[1][:st0_stop]**n/(1+hd*un_dc[1][:st0_stop]**n),color[2]+'-',label='death feedback, no treatment')
    axD.legend();
    axD.set_ylabel("Death Feedback (relative)")
    axD.set_title(s_title)
    if saveQ:
        figD.savefig(plot_drty+"\\death_fdbk_"+plot_save_suffix+".png",dpi=300)

    figD, axD = plt.subplots(figsize=(8,8));
    if log_yscaleQ:
        axD.set_yscale('log')
    axD.plot(t_vec[0][:st2_stop],1/(1+l*u_dc[0][:st2_stop]**n),color[0]+'-',label='prob feedback '+str(csr[0])+unit)
    axD.plot(t_vec[1][:st1_stop],1/(1+l*u_dc[1][:st1_stop]**n),color[1]+'-',label='prob feedback '+str(csr[1])+unit)
    axD.plot(tn_vec[0][:st0_stop],1/(1+l*un_dc[1][:st0_stop]**n),color[2]+'-',label='prob feedback, no treatment')
    ### fixed to actually being feedback on renewal and not a copy-paste of feedback on death as of 30 June...so all plots before this date are copies of the feedback on death one
    axD.legend();
    axD.set_ylabel("Probability Feedback (relative)")
    axD.set_title(s_title)
    if saveQ:
        figD.savefig(plot_drty+"\\prob_fdbk_"+plot_save_suffix+".png",dpi=300)

    figR1, axR1 = plt.subplots(figsize=(8,8));
    if log_yscaleQ:
        axR1.set_yscale('log')
    axR1.plot(t_vec[0][:st2_stop],1/(1+h1*u_dc[0][:st2_stop]**n),color[0]+'-',label='csc division rate feedback '+str(csr[0])+unit)
    axR1.plot(t_vec[1][:st1_stop],1/(1+h1*u_dc[1][:st1_stop]**n),color[1]+'-',label='csc division rate feedback '+str(csr[1])+unit)
    axR1.plot(tn_vec[0][:st0_stop],1/(1+h1*un_dc[1][:st0_stop]**n),color[2]+'-',label='csc division rate feedback, no treatment')
    axR1.legend();
    axR1.set_ylabel("CSC Division Rate Feedback (relative)")
    axR1.set_title(s_title)
    if saveQ:
        figR1.savefig(plot_drty+"\\csc_div_fdbk_"+plot_save_suffix+".png",dpi=300)

    figR2, axR2 = plt.subplots(figsize=(8,8));
    if log_yscaleQ:
        axR2.set_yscale('log')
    axR2.plot(t_vec[0][:st2_stop],1/(1+h2*u_dc[0][:st2_stop]**n),color[0]+'-',label='dcc division rate feedback '+str(csr[0])+unit)
    axR2.plot(t_vec[1][:st1_stop],1/(1+h2*u_dc[1][:st1_stop]**n),color[1]+'-',label='dcc division rate feedback '+str(csr[1])+unit)
    axR2.plot(tn_vec[0][:st0_stop],1/(1+h2*un_dc[1][:st0_stop]**n),color[2]+'-',label='dcc division rate feedback, no treatment')
    axR2.legend();
    axR2.set_ylabel("DCC Division Rate Feedback (relative)")
    axR2.set_title(s_title)
    if saveQ:
        figR2.savefig(plot_drty+"\\dcc_div_fdbk_"+plot_save_suffix+".png",dpi=300)
        
    figEffR2, axEffR2 = plt.subplots(figsize=(8,8));
    if log_yscaleQ:
        axR2.set_yscale('log')
    axEffR2.plot(t_vec[0][:st2_stop],r2/(1+h2*u_dc[0][:st2_stop]**n) - r2/5-(1.1*r2-r2/5)*hd*u_dc[0][:st2_stop]**n/(1+hd*u_dc[0][:st2_stop]**n),color[0]+'-',label='effective dcc division rate '+str(csr[0])+unit)
    axEffR2.plot(t_vec[1][:st1_stop],r2/(1+h2*u_dc[1][:st1_stop]**n) - r2/5-(1.1*r2-r2/5)*hd*u_dc[1][:st1_stop]**n/(1+hd*u_dc[1][:st1_stop]**n),color[1]+'-',label='ffective dcc division rate '+str(csr[1])+unit)
    axEffR2.plot(tn_vec[0][:st0_stop],r2/(1+h2*un_dc[1][:st0_stop]**n) - r2/5-(1.1*r2-r2/5)*hd*un_dc[1][:st0_stop]**n/(1+hd*un_dc[1][:st0_stop]**n),color[2]+'-',label='ffective dcc division rate, no treatment')
    axEffR2.legend();
    axEffR2.set_ylabel("Effective DCC Division Rate")
    axEffR2.set_title(s_title)
    if saveQ:
        figEffR2.savefig(plot_drty+"\\dcc_eff_div_fdbk_"+plot_save_suffix+".png",dpi=300)
        
    n_t = -180000*2; #-2200000; #-180000
    figCa, axCa = plt.subplots(figsize=(8,8));
    if log_yscaleQ:
        axCa.set_yscale('log')
    axCa.plot(t_vec[0][n_t:],csc_p[0][n_t:]/csc_p[1][n_t:],color[3]+':',label='CSC '+str(csr[0])+':' +str(csr[1]) + unit)
    axCa.legend();
    axCa.set_ylabel("CSC frac ratio")
    axCa.set_title('t_vec[0]['+str(n_t)+':];'+s_title)
    if saveQ:
        print("saving graph of index",lll_load)
        figCa.savefig(plot_drty+"\\csc_ratio_"+plot_save_suffix+".png",dpi=300)


# TODO: Check what this produces against what Nayeong has produced.