# -*- coding: utf-8 -*-
"""
Created on Sun May 23 12:50:58 2021

@author: jhvo9
"""
import matplotlib.pyplot as plt;
import numpy as np;
total_cell_num = 4/3*np.pi*(1 ** 3)*1e9;
drty = "C:\\Users\\jhvo9\\Google Drive (vojh1@uci.edu)\\a PhD Projects\\GBM Modeling\\python scripts\\data\\kim_model\\22_May\\"
TU_2 = np.loadtxt(drty+"\\with_death_fdbk\\U_2Gy.txt");
TU_1 = np.loadtxt(drty+"\\with_death_fdbk\\U_1Gy.txt");
TU_2_none = np.loadtxt(drty+"\\with_death_fdbk\\U_none_2Gy.txt");
TU_1_none = np.loadtxt(drty+"\\with_death_fdbk\\U_none_1Gy.txt");
t_vec = [TU_2[0,:], TU_1[0,:]];
u_sc = [TU_2[1,:], TU_1[1,:]];
u_dc = [TU_2[2,:], TU_1[2,:]];
u_srv = [TU_2[3,:], TU_1[3,:]];

tn_vec = [TU_2_none[0,:], TU_1_none[0,:]];
un_sc = [TU_2_none[1,:], TU_1_none[1,:]];
un_dc = [TU_2_none[2,:], TU_1_none[2,:]];
un_srv = [TU_2_none[3,:], TU_1_none[3,:]];
Doses = [2,1];
hd = 1e5; n = 1; deathFdbkQ = True;
figT, axT = plt.subplots(figsize=(8,8));
axT.set_yscale('log')
axT.plot(t_vec[0],u_sc[0] * total_cell_num,'bo--',label='CSC '+str(Doses[0])+' Gy')
axT.plot(t_vec[0],u_dc[0] * total_cell_num,'b-',label='DCC '+str(Doses[0])+' Gy')
axT.plot(t_vec[1],u_sc[1] * total_cell_num,'go--',label='CSC '+str(Doses[1])+' Gy')
axT.plot(t_vec[1],u_dc[1] * total_cell_num,'g-',label='DCC '+str(Doses[1])+' Gy')
axT.legend();
axT.set_ylabel("Cell Number (log)")
plt.savefig(figT,"")

# TODO: maybe add the other ones too?
csc_p = [];
csc_p += [u_sc[0]/(u_sc[0] + u_dc[0])];
csc_p += [u_sc[1]/(u_sc[1] + u_dc[1])];
figC, axC = plt.subplots(figsize=(8,8));
axC.set_yscale('log')
axC.plot(t_vec[0],csc_p[0],'bo--',label='CSC '+str(Doses[0])+' Gy')
axC.plot(t_vec[1],csc_p[1],'go--',label='CSC '+str(Doses[1])+' Gy')
axC.legend();
axC.set_ylabel("CSC %")

figD, axD = plt.subplots(figsize=(8,8));
axD.set_yscale('log')
axD.plot(t_vec[0],1/(1+hd*u_dc[0]**n),'bo--',label='death feedback '+str(Doses[0])+' Gy')
axD.plot(t_vec[1],1/(1+hd*u_dc[1]**n),'go--',label='death feedback '+str(Doses[1])+' Gy')
axD.legend();
axD.set_ylabel("Death Feedback (relative)")

n = -2200000;
figCa, axCa = plt.subplots(figsize=(8,8));
axCa.set_yscale('log')
axCa.plot(t_vec[0][n:],csc_p[0][n:]/csc_p[1][n:],'r--',label='CSC '+str(Doses[0])+':' +str(Doses[1]) + 'Gy')
axCa.legend();
axCa.set_ylabel("CSC % ratio (2 Gy to 1 Gy)")
axCa.set_title('t_vec[0]['+str(n)+':]; Death Feedback:'+str(deathFdbkQ))
# TODO: Check what this produces against what Nayeong has produced.