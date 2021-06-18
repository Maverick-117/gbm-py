# -*- coding: utf-8 -*-
"""
Created on Mon May 24 14:20:07 2021

@author: jhvo9
"""
import matplotlib.pyplot as plt;
import numpy as np;
total_cell_num = 4/3*np.pi*(1 ** 3)*1e9;
compDosesQ = False; #true: yes you're comparing doses; false: no, you're comparing reprogramming versus non-reprogramming
saveQ = True;
Doses = [2,1];
C = [0, 0.005196];
deathFdbkQ = True;
if compDosesQ:
    unit = ' Gy';
    csr = Doses;
else:
    unit = ' reprog'
    csr = C;
lll_load = 4;
color = ['k','r','g','b','m'];
if compDosesQ:
    ending1 = str(Doses[0])+unit[1:];
    ending2 = str(Doses[1])+unit[1:];
else:
    ending1 = str(C[0])+unit[1:];
    ending2 = str(C[1])+unit[1:];
if deathFdbkQ:
    t1 = "_w_death_fdbk";
    ending1+=t1;
    ending2+=t1;
else:
    t2 = "_w_no_death_fdbk";
    ending1+=t2;
    ending2+=t2;
    
plot_str2, plot_str1 = "TU_"+str(csr[0])+unit[1:]+"_"+str(lll_load), "TU_"+str(csr[1])+unit[1:]+"_"+str(lll_load);
plotn_str2, plotn_str1 = "TU_none_"+str(csr[0])+unit[1:]+"_"+str(lll_load), "TU_none_"+str(csr[1])+unit[1:]+"_"+str(lll_load);

data_drty = "C:\\Users\\jhvo9\\Google Drive (vojh1@uci.edu)\\a PhD Projects\\GBM Modeling\\python scripts\\data\\kim_model\\25_May\\kim_reprog\\"
plot_drty = "C:\\Users\\jhvo9\\Google Drive (vojh1@uci.edu)\\a PhD Projects\\GBM Modeling\\python scripts\\plots\\kim_model\\25_May\\kim_reprog\\"
TU_2 = np.loadtxt(data_drty+plot_str2+".txt");
TU_1 = np.loadtxt(data_drty+plot_str1+".txt");
TU_2_none = np.loadtxt(data_drty+plotn_str2+".txt");
TU_1_none = np.loadtxt(data_drty+plotn_str1+".txt");

hd = 1e5; n = 1; 

t_vec = [TU_2[0,:], TU_1[0,:]];
u_sc = [TU_2[1,:], TU_1[1,:]];
u_dc = [TU_2[2,:], TU_1[2,:]];
u_srv = [TU_2[3,:], TU_1[3,:]];

tn_vec = [TU_2_none[0,:], TU_1_none[0,:]];
un_sc = [TU_2_none[1,:], TU_1_none[1,:]];
un_dc = [TU_2_none[2,:], TU_1_none[2,:]];
un_srv = [TU_2_none[3,:], TU_1_none[3,:]];

T_stop = 200;
#T_stop = int(t_vec[0].max());

st2_stop = np.flatnonzero(t_vec[0] < T_stop).max();
st1_stop = np.flatnonzero(t_vec[1] < T_stop).max();
st0_stop = np.flatnonzero(tn_vec[0] < T_stop).max();

plot_save_suffix = '_'+str(T_stop)+'days_'+str(lll_load);

Doses = [2,1];
C = [0, 0.005196];
hd = 1e5; n = 1; deathFdbkQ = True;
figT, axT = plt.subplots(figsize=(8,8));
axT.set_yscale('log')
axT.plot(t_vec[0][:st2_stop],u_sc[0][:st2_stop] * total_cell_num,color[0]+'--',label='CSC '+str(csr[0])+unit)
axT.plot(t_vec[0][:st2_stop],u_dc[0][:st2_stop] * total_cell_num,color[0]+'-',label='DCC '+str(csr[0])+unit)
axT.plot(t_vec[1][:st1_stop],u_sc[1][:st1_stop] * total_cell_num,color[1]+'--',label='CSC '+str(csr[1])+unit)
axT.plot(t_vec[1][:st1_stop],u_dc[1][:st1_stop] * total_cell_num,color[1]+'-',label='DCC '+str(csr[1])+unit)
axT.plot(tn_vec[0][:st0_stop],un_sc[0][:st0_stop] * total_cell_num,color[2]+'--',label='CSC, no treatment')
axT.plot(tn_vec[0][:st0_stop],un_dc[0][:st0_stop] * total_cell_num,color[2]+'-',label='DCC, no treatment')
axT.legend();
axT.set_ylabel("Cell Number (log)")
axT.set_title("Death Feedback:"+str(deathFdbkQ))
if saveQ:
    figT.savefig(plot_drty+"\\pop_plot_"+plot_save_suffix+".png",dpi=300)

# TODO: maybe add the other ones too?
csc_p = [];
csc_p += [u_sc[0]/(u_sc[0] + u_dc[0])];
csc_p += [u_sc[1]/(u_sc[1] + u_dc[1])];
csc_p += [un_sc[0]/(un_sc[0] + un_dc[0])];
figC, axC = plt.subplots(figsize=(8,8));
axC.set_yscale('log')
axC.plot(t_vec[0][:st2_stop],csc_p[0][:st2_stop],color[0]+'-',label='CSC '+str(csr[0])+unit)
axC.plot(t_vec[1][:st1_stop],csc_p[1][:st1_stop],color[1]+'-',label='CSC '+str(csr[1])+unit)
axC.plot(tn_vec[0][:st0_stop],csc_p[2][:st0_stop],color[2]+'-',label='CSC, no treatment')
axC.legend();
axC.set_ylabel("CSC %")
axC.set_title("Death Feedback:"+str(deathFdbkQ))
if saveQ:
    figC.savefig(plot_drty+"\\csc_plot_"+plot_save_suffix+".png",dpi=300)

figD, axD = plt.subplots(figsize=(8,8));
axD.set_yscale('log')
axD.plot(t_vec[0][:st2_stop],1/(1+hd*u_dc[0][:st2_stop]**n),color[0]+'-',label='death feedback '+str(csr[0])+unit)
axD.plot(t_vec[1][:st1_stop],1/(1+hd*u_dc[1][:st1_stop]**n),color[1]+'-',label='death feedback '+str(csr[1])+unit)
axD.plot(tn_vec[0][:st0_stop],1/(1+hd*un_dc[1][:st0_stop]**n),color[2]+'-',label='death feedback, no treatment')
axD.legend();
axD.set_ylabel("Death Feedback (relative)")
axD.set_title("Death Feedback:"+str(deathFdbkQ))
if saveQ:
    figD.savefig(plot_drty+"\\death_fdbk_"+plot_save_suffix+".png",dpi=300)

n = -180000; #-2200000; #-180000
figCa, axCa = plt.subplots(figsize=(8,8));
axCa.set_yscale('log')
axCa.plot(t_vec[0][n:],csc_p[0][n:]/csc_p[1][n:],color[3]+':',label='CSC '+str(csr[0])+':' +str(csr[1]) + unit)
axCa.legend();
axCa.set_ylabel("CSC % ratio (2 Gy to 1 Gy)")
axCa.set_title('t_vec[0]['+str(n)+':]; Death Feedback:'+str(deathFdbkQ))
if saveQ:
    figCa.savefig(plot_drty+"\\csc_ratio_"+plot_save_suffix+".png",dpi=300)

# TODO: Check what this produces against what Nayeong has produced.