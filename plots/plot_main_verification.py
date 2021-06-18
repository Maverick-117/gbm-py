# -*- coding: utf-8 -*-
"""
Created on Mon May 24 14:20:07 2021

@author: jhvo9
"""
import matplotlib.pyplot as plt;
import numpy as np;
total_cell_num = 4/3*np.pi*(1 ** 3)*1e9;
compDosesQ = False; #true: yes you're comparing doses; false: no, you're comparing reprogramming versus non-reprogramming
saveQ = False;
Doses = [2,1];
C = [0, 0.005196];
if compDosesQ:
    unit = ' Gy';
    csr = Doses;
else:
    unit = ' reprog'
    csr = C;
lll_load = 0;
plot_str2, plot_str1 = "TU_"+str(csr[0])+unit[1:]+"_"+str(lll_load), "TU_"+str(csr[1])+unit[1:]+"_"+str(lll_load);
plotn_str2, plotn_str1 = "TU_none_"+str(csr[0])+unit[1:]+"_"+str(lll_load), "TU_none_"+str(csr[1])+unit[1:]+"_"+str(lll_load);

data_drty = "C:\\Users\\jhvo9\\Google Drive (vojh1@uci.edu)\\a PhD Projects\\GBM Modeling\\python scripts\\data\\kim_model\\verification\\"
plot_drty = "C:\\Users\\jhvo9\\Google Drive (vojh1@uci.edu)\\a PhD Projects\\GBM Modeling\\python scripts\\plots\\kim_model\\verification\\"
TU_2 = np.loadtxt(data_drty+plot_str2+".txt");
TU_1 = np.loadtxt(data_drty+plot_str1+".txt");
TU_2_none = np.loadtxt(data_drty+plotn_str2+".txt");
TU_1_none = np.loadtxt(data_drty+plotn_str1+".txt");

hd = 1e5; n = 1; deathFdbkQ = True;

t_vec = [TU_2[0,:], TU_1[0,:]];
u_sc = [TU_2[1,:], TU_1[1,:]];
u_dc = [TU_2[2,:], TU_1[2,:]];
u_srv = [TU_2[3,:], TU_1[3,:]];

tn_vec = [TU_2_none[0,:], TU_1_none[0,:]];
un_sc = [TU_2_none[1,:], TU_1_none[1,:]];
un_dc = [TU_2_none[2,:], TU_1_none[2,:]];
un_srv = [TU_2_none[3,:], TU_1_none[3,:]];

st2_500 = np.flatnonzero(t_vec[0] < 500).max();
st1_500 = np.flatnonzero(t_vec[1] < 500).max();

Doses = [2,1];
C = [0, 0.005196];
hd = 1e5; n = 1; deathFdbkQ = True;
figT, axT = plt.subplots(figsize=(8,8));
axT.set_yscale('log')
axT.plot(t_vec[0],u_sc[0] * total_cell_num,'bo--',label='CSC '+str(csr[0])+unit)
axT.plot(t_vec[0],u_dc[0] * total_cell_num,'b-',label='DCC '+str(csr[0])+unit)
axT.plot(t_vec[1],u_sc[1] * total_cell_num,'go--',label='CSC '+str(csr[1])+unit)
axT.plot(t_vec[1],u_dc[1] * total_cell_num,'g-',label='DCC '+str(csr[1])+unit)
axT.legend();
axT.set_ylabel("Cell Number (log)")
if saveQ:
    figT.savefig(plot_drty+"\\pop_plot"+str(lll_load)+".png",dpi=300)

# TODO: maybe add the other ones too?
csc_p = [];
csc_p += [u_sc[0]/(u_sc[0] + u_dc[0])];
csc_p += [u_sc[1]/(u_sc[1] + u_dc[1])];
figC, axC = plt.subplots(figsize=(8,8));
axC.set_yscale('log')
axC.plot(t_vec[0][:st2_500],csc_p[0][:st2_500],'bo--',label='CSC '+str(csr[0])+unit)
axC.plot(t_vec[1][:st1_500],csc_p[1][:st1_500],'go--',label='CSC '+str(csr[1])+unit)
axC.legend();
axC.set_ylabel("CSC %")

figD, axD = plt.subplots(figsize=(8,8));
axD.set_yscale('log')
axD.plot(t_vec[0],1/(1+hd*u_dc[0]**n),'bo--',label='death feedback '+str(csr[0])+unit)
axD.plot(t_vec[1],1/(1+hd*u_dc[1]**n),'go--',label='death feedback '+str(csr[1])+unit)
axD.legend();
axD.set_ylabel("Death Feedback (relative)")

n = -2200000; #-180000
figCa, axCa = plt.subplots(figsize=(8,8));
axCa.set_yscale('log')
axCa.plot(t_vec[0][n:],csc_p[0][n:]/csc_p[1][n:],'r--',label='CSC '+str(csr[0])+':' +str(csr[1]) + unit)
axCa.legend();
axCa.set_ylabel("CSC % ratio (2 Gy to 1 Gy)")
axCa.set_title('t_vec[0]['+str(n)+':]; Death Feedback:'+str(deathFdbkQ))
# TODO: Check what this produces against what Nayeong has produced.