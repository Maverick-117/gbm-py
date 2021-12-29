# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 23:32:13 2021

Fix reprogramming level. Vary feedback types over the same graph.
now features a vary-over-mu-relaxation-rate loop
3 Nov: now features a vary-over-beta loop
22 Nov: now does variation over schedule type

Created on Wed Nov  3 20:41:24 2021

@author: jhvo9
"""
import matplotlib.pyplot as plt;
import numpy as np;
from os import makedirs, chdir, replace
from os.path import exists
import matplotlib;
from funciones import calc_EOT, calc_BED_Frac, dU_dt;
from numpy import sqrt
import matplotlib as mpl
from matplotlib import cm

class MplColorHelper:

  def __init__(self, cmap_name, start_val, stop_val):
    self.cmap_name = cmap_name
    self.cmap = plt.get_cmap(cmap_name)
    self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
    self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

  def get_rgb(self, val):
    return self.scalarMap.to_rgba(val)
def calc_eigs6(U,V):
    #assumes: using feedback type 6, with feedback on death, k2 model, l = 0, h1 = 1e5, h2 = 0, hd = h1/3, c = 5.169e-3
    eig1 = (1/(2*(1. + 100000.*V)**2*(123. + 4.e6*V)**2))*((-(1/(1. + 100000.*V)))*(-26.88877870433699 -7.126619396434031e6*V - 6.47097325564284e11*V**2 -2.3175998037183836e16*V**3 - 2.8436807407587528e20*V**4 -2151.1022963469572/(1. + 100000.*V) +(2.66198909172936e8*U)/(1. +100000.*V) - (8.429522819831166e8*V)/(1. + 100000.*V) +(7.055353202469684e13*U*V)/(1. +100000.*V) - (1.2178631508447506e14*V**2)/(1. + 100000.*V) +(6.406263523086407e18*U*V**2)/(1. +100000.*V) - (7.441343762417502e18*V**3)/(1. + 100000.*V) +(2.2944238056811978e23*U*V**3)/(1. +100000.*V) - (1.2809359896747807e23*V**4)/(1. + 100000.*V) +(2.815243933351163e27*U*V**4)/(1. +100000.*V) + (2.886335951870129e27*V**5)/(1. + 100000.*V) +(2.8436807407587454e31*V**6)/(1. + 100000.*V)) -sqrt((-(1/(1. + 100000.*V)**3))*4*(57840.51361686446 +3.799604471741991e10*V - 0.5*U*V +1.0674057901553622e16*V**2 +1.663339470839293e21*V**3 + 1.5543146926018233e26*V**4 +8.712766041868346e30*V**5 +2.686472848142641e35*V**6 + 3.014379014122504e39*V**7 -4.886936833640244e43*V**8 -1.6225927682921336e32*U*V**8 - 1.4798331884312886e48*V**9 -8.086520155362228e51*V**10) + (1/(1. + 100000.*V)**2)*(-26.88877870433699 -7.126619396434031e6*V -6.47097325564284e11*V**2 -2.3175998037183836e16*V**3 - 2.8436807407587528e20*V**4 -2151.1022963469572/(1. + 100000.*V) + (2.66198909172936e8*U)/(1. + 100000.*V) -(8.429522819831166e8*V)/(1. +100000.*V) + (7.055353202469684e13*U*V)/(1. +100000.*V) -(1.2178631508447506e14*V**2)/(1. +100000.*V) + (6.406263523086407e18*U*V**2)/(1. +100000.*V) -(7.441343762417502e18*V**3)/(1. +100000.*V) + (2.2944238056811978e23*U*V**3)/(1. +100000.*V) -(1.2809359896747807e23*V**4)/(1. +100000.*V) + (2.815243933351163e27*U*V**4)/(1. +100000.*V) +(2.886335951870129e27*V**5)/(1. +100000.*V) + (2.8436807407587454e31*V**6)/(1. +100000.*V))**2))
    eig2 = (1/(2*(1. + 100000.*V)**2*(123. + 4.e6*V)**2))*((-(1/(1. + 100000.*V)))*(-26.88877870433699 -7.126619396434031e6*V - 6.47097325564284e11*V**2 -2.3175998037183836e16*V**3 - 2.8436807407587528e20*V**4 -2151.1022963469572/(1. + 100000.*V) +(2.66198909172936e8*U)/(1. +100000.*V) - (8.429522819831166e8*V)/(1. + 100000.*V) +(7.055353202469684e13*U*V)/(1. +100000.*V) - (1.2178631508447506e14*V**2)/(1. + 100000.*V) +(6.406263523086407e18*U*V**2)/(1. +100000.*V) - (7.441343762417502e18*V**3)/(1. + 100000.*V) +(2.2944238056811978e23*U*V**3)/(1. +100000.*V) - (1.2809359896747807e23*V**4)/(1. + 100000.*V) +(2.815243933351163e27*U*V**4)/(1. +100000.*V) + (2.886335951870129e27*V**5)/(1. + 100000.*V) +(2.8436807407587454e31*V**6)/(1. + 100000.*V)) +sqrt((-(1/(1. + 100000.*V)**3))*4*(57840.51361686446 +3.799604471741991e10*V - 0.5*U*V +1.0674057901553622e16*V**2 +1.663339470839293e21*V**3 + 1.5543146926018233e26*V**4 +8.712766041868346e30*V**5 +2.686472848142641e35*V**6 + 3.014379014122504e39*V**7 -4.886936833640244e43*V**8 -1.6225927682921336e32*U*V**8 - 1.4798331884312886e48*V**9 -8.086520155362228e51*V**10) + (1/(1. + 100000.*V)**2)*(-26.88877870433699 -7.126619396434031e6*V -6.47097325564284e11*V**2 -2.3175998037183836e16*V**3 - 2.8436807407587528e20*V**4 -2151.1022963469572/(1. + 100000.*V) + (2.66198909172936e8*U)/(1. + 100000.*V) -(8.429522819831166e8*V)/(1. +100000.*V) + (7.055353202469684e13*U*V)/(1. +100000.*V) -(1.2178631508447506e14*V**2)/(1. +100000.*V) + (6.406263523086407e18*U*V**2)/(1. +100000.*V) -(7.441343762417502e18*V**3)/(1. +100000.*V) + (2.2944238056811978e23*U*V**3)/(1. +100000.*V) -(1.2809359896747807e23*V**4)/(1. +100000.*V) + (2.815243933351163e27*U*V**4)/(1. +100000.*V) +(2.886335951870129e27*V**5)/(1. +100000.*V) + (2.8436807407587454e31*V**6)/(1. +100000.*V))**2))
    return [eig1, eig2];

font = {'weight' : 'bold',
        'size'   : 13}

l_w = 10**(-7); # weak feedback on prob
l_s = 10**3; # strong feedback on prob
h1 = 10**5; # feedback on css div 
h2 = 10**5; # feedback on dcc div
l_vec =  [0,  0, l_w, l_s, l_w , l_s , 0 ,0 ,l_w,l_w,l_s,l_s];
h1_vec = [0, h1, 0  , 0  , h1  , h1  , h1,0 ,h1 ,0  , h1, 0 ];
h2_vec = [0, h2, 0  , 0  , h2  , h2  , 0 ,h2,0  ,h2 , 0 , h2];


load_temp = np.arange(12).tolist();#[1,3,5,6,7,10,11];#np.arange(12).tolist();


#[1,6,7] had...the same stability behaviour/eigenvalue plots
saveQ = False;
cases = ["2 Gy", "2.67 Gy" , "3.4 Gy", "5 Gy"]; 
comp_str = "Gy";
Doses = [2.0,6.0];
Fracs = [calc_BED_Frac(.17,.02,Doses[0]),calc_BED_Frac(.17,.02,Doses[1])]#list(map(lambda d:calc_BED_Frac(a,b,d,40*(1+(3.4)/8.5)),Doses));#;[25,20]
EOTs  = [calc_EOT(Fracs[0],100), calc_EOT(Fracs[1],100)];
EOT_diff = EOTs[0] - EOTs[1];
EOT_str = ",EOT-shifted"
match_EOT =True;
if not match_EOT:
    EOT_diff = 0.;
    EOT_str = "";
C = [5.196*10**(-3), 5.196*10**(-3)];
comp_dir = str(C[0])+"_reprog";
fix_str= str(Doses[-1])+"Gy"
comp_list = C;
unit = ' Gy';
T_stop = 1200;
#lll_load = 6;

deathFdbkQ = False;
if deathFdbkQ:
    hd = 32520.32520325203;#1e5; 
    deathFdbk_str = "_w_death_fdbk";
else:
    hd = 0.0;#1e5; 
    deathFdbk_str = "_w_no_death_fdbk";
sub_drty = "\\corrected_reprog";
deathVal_dir = '\\death_val_of_kim';
color = ['r','k','g','b','m'];


compBaselineQ = False;
showNonReprogQ = True;
schedChangeQ = False;
use_muQ = True;
if use_muQ:
    #, 1e3,np.infty]
    # sig_vec = [10];#list(map(lambda p: 1 * 10 ** (p), np.arange(0,3).tolist()));
    # rho_vec = [0];#[1/5];#sorted(list(map(lambda p: 2 ** (p), np.arange(0,-5,-1).tolist())) + [1/5]);
    sig_vec = [0.1];#sorted(list(map(lambda p:  10 ** (p), np.arange(-1,3).tolist()))+[0]);
    rho_vec = [0.2];#sorted(list(map(lambda p: 2 ** (p), np.arange(0,-3,-1).tolist())) + [0,2]);
# rename rho to rho
    xi1_vec = [0.01];#sorted(list(map(lambda p:  10 ** (p), np.arange(-1,3).tolist()))+[0]);
    xi2_vec = [1];#sorted(list(map(lambda p:  10 ** (p), np.arange(-1,3).tolist()))+[0]);#np.arange(0.01,0.1,0.01);#[.1];#sorted(list(map(lambda p:  10 ** (p), np.arange(0,3).tolist()))+[0])
    s_xi1_vec = ['\\' + str(el) for el in xi1_vec];
    s_xi2_vec = ['\\' + str(elb) for elb in xi2_vec];
    s_sig_vec = ['\\' + str(el) for el in sig_vec];
    s_rho_vec = ['\\' + str(elb) + '\\' for elb in rho_vec];
    case = 'feedbackCheck\\21_Dec'#'reversionAttempt4\\14_Nov';#'muNN-l2b\\2_Nov';#'baseline\\31_Oct';
else:
    sig_vec=[0];
    s_sig_vec = ['\\0\\'];
    case = 'muNN\\26_Oct';#'full_fdbk_gains\\13_Oct'

if schedChangeQ:
    Nsched_vec = [1];
else:
    Nsched_vec = [0];
# Nsched_vec = [0,1];
Nsched_vec_str = ['',' (delayed)']
base_dir = "G:\\My Drive\\a PhD Projects\\GBM Modeling\\python scripts\\data\\k2_model\\"
base_dir1 = base_dir+"conventional_BED\\"+case+"\\death_val_of_kim\\corrected_reprog";#"conventional_BED\\srv_test\\16_Oct\\death_val_of_kim\\corrected_reprog";
base_dir2 = base_dir+"fixed_BED\\40\\03_Aug\\death_val_of_kim\\corrected_reprog"




data2 = np.ones((len(Nsched_vec),len(xi1_vec),len(xi2_vec),len(rho_vec),len(sig_vec),12)).tolist();
data1 = np.ones((len(Nsched_vec),len(xi1_vec),len(xi2_vec),len(rho_vec),len(sig_vec),12)).tolist();
dataN = np.ones((len(Nsched_vec),len(xi1_vec),len(xi2_vec),len(rho_vec),len(sig_vec),12)).tolist();
for lll_load in load_temp:
    for nn,eln in enumerate(Nsched_vec):
        s_fr_0 = '_'+str(Fracs[0])+'_Fracs_';
        s_fr_1 = '_'+str(Fracs[1])+'_Fracs_';#'_'+str(Fracs[1])+'_Fracs_';

        for tt, elt in enumerate(s_xi1_vec):
            for uu, elu in enumerate(s_xi2_vec):
                for ww, elb in enumerate(s_rho_vec):
                    for vv, el in enumerate(s_sig_vec):
                        s_mu_str = elt + elu  + elb  + el +"\\Schedule"+str(eln) + "\\";
                        plot_str2 = "TU_"+str(C[0])+'_reprog'+str(Doses[0])+'_'+comp_str+s_fr_0+str(T_stop)+'days_'+deathFdbk_str+"_"+str(lll_load);
                        plot_str1 = "TU_"+str(C[1])+'_reprog'+str(Doses[1])+'_'+comp_str+s_fr_1+str(T_stop)+'days_'+deathFdbk_str+"_"+str(lll_load);
                        data2[nn][tt][uu][ww][vv][lll_load] = np.loadtxt(base_dir1 + s_mu_str + plot_str2 + ".txt");
                        data1[nn][tt][uu][ww][vv][lll_load] = np.loadtxt(base_dir1 + s_mu_str + plot_str1 + ".txt");
                        print('dir2',base_dir1 + s_mu_str + plot_str2 + ".txt")
                        print('dir1',base_dir1 + s_mu_str + plot_str1 + ".txt")
                        plotn_str2 = "TU_none_"+str(C[1])+'_reprog'+str(Doses[0])+'_'+comp_str+s_fr_0+str(T_stop)+'days_'+deathFdbk_str+"_"+str(lll_load); 
                        print('dir_none',base_dir1 + s_mu_str + plotn_str2 + ".txt")
                        dataN[nn][tt][uu][ww][vv][lll_load] = np.loadtxt(base_dir1 + s_mu_str + plotn_str2 + ".txt");
log_yscaleQ = True;
if log_yscaleQ:
    log_str = " (log)";
    log_dir = "\\log_yscale\\";
else:
    log_str = "";
    log_dir = "\\linear_yscale\\"
#T_stop = 300;
T_stop = 1100;
total_cell_num = 4/3*np.pi*(1 ** 3)*1e9;
s_sched0 = str((Doses[0],Fracs[0]));
s_sched1 = str((Doses[1],Fracs[1]));
s_reprg0 = f"{float(C[0]):.3}";
s_reprg1 = f"{float(C[1]):.3}";
s_reprg_title = ";c:";


if use_muQ:
    s_sched_title = ";sched:"
    s_title_tweakable = s_sched_title+s_sched0;
    s_lgd_0 = s_sched0;#str(csr[0])+unit+EOT_str;
    s_lgd_1 = s_sched1+EOT_str;#str(csr[1])+unit;
    COL2 = MplColorHelper('viridis', 0, len(Nsched_vec)-1)
    COL1 = MplColorHelper('cool', 0, len(Nsched_vec)-1)
    for lll_load in load_temp:
        figTcsc, axTcsc = plt.subplots(figsize=(8,8));
        figTdcc, axTdcc = plt.subplots(figsize=(8,8));
        s_title = "(l:"+f"{float(l_vec[lll_load]):.3}"+";h1:"+f"{float(h1_vec[lll_load]):.3}"+";h2:"+f"{float(h2_vec[lll_load]):.3}"+";hd:"+f"{hd:.3}"+s_title_tweakable+")";
        # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
        # s_lgd_1 = "w/ constant reprog";#"Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';

        for nn,eln in enumerate(Nsched_vec):
            # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
            # s_lgd_1 = "Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
            for tt, elt in enumerate(xi1_vec):
                xi1 = elt;
                s_xi1 = '('+str(elt)+',';
                for uu, elu in enumerate(xi2_vec):
                    xi2 = elu;
                    s_xi2 = str(elu)+',';
                    for ww, elb in enumerate(rho_vec):
                        beta = elb;
                        s_beta = str(elb) + ',';
                        for vv, el in enumerate(sig_vec):
                            s_sig = '(σ='+str(el) +')';
                            s_title_tweakable = s_xi1+s_xi2+s_beta+s_sig;
                            st2_stop = np.flatnonzero(data2[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
                            axTcsc.plot(data2[nn][tt][uu][ww][vv][lll_load][0,:st2_stop],data2[nn][tt][uu][ww][vv][lll_load][1,:st2_stop] * total_cell_num,'--',color=COL2.get_rgb(nn),label='CSC '+s_lgd_0+s_sig+'_'+str(lll_load))
                            axTdcc.plot(data2[nn][tt][uu][ww][vv][lll_load][0,:st2_stop],data2[nn][tt][uu][ww][vv][lll_load][2,:st2_stop] * total_cell_num,'-',color=COL2.get_rgb(nn),label='DCC '+s_lgd_0+s_sig+'_'+str(lll_load))
                            if showNonReprogQ:
                                st1_stop = np.flatnonzero(data1[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
                                axTcsc.plot(data1[nn][tt][uu][ww][vv][lll_load][0,:st2_stop]+EOT_diff,data1[nn][tt][uu][ww][vv][lll_load][1,:st2_stop] * total_cell_num,'--',color=COL1.get_rgb(nn),label='CSC '+s_lgd_1+s_sig+'_'+str(lll_load))
                                axTdcc.plot(data1[nn][tt][uu][ww][vv][lll_load][0,:st1_stop]+EOT_diff,data1[nn][tt][uu][ww][vv][lll_load][2,:st1_stop] * total_cell_num,'-',color=COL1.get_rgb(nn),label='DCC '+s_lgd_1+s_sig+'_'+str(lll_load))
                        
        st0_stop = np.flatnonzero(dataN[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
        axTcsc.plot(dataN[nn][tt][uu][ww][vv][lll_load][0,:st0_stop],dataN[nn][tt][uu][ww][vv][lll_load][1,:st0_stop] * total_cell_num,'g--',label='CSC, no treatment_'+str(lll_load))
        axTdcc.plot(dataN[nn][tt][uu][ww][vv][lll_load][0,:st0_stop],dataN[nn][tt][uu][ww][vv][lll_load][2,:st0_stop] * total_cell_num,'g-',label='DCC, no treatment_'+str(lll_load))
        axTcsc.legend(fancybox=True, framealpha=0.5);
        axTcsc.set_ylabel("CSC Cell Number"+log_str)
        axTcsc.set_xlabel("time (days)")
        axTcsc.set_title(s_title)
        axTdcc.legend(fancybox=True, framealpha=0.5);
        axTdcc.set_ylabel("DCC Cell Number"+log_str)
        axTdcc.set_xlabel("time (days)")
        axTdcc.set_title(s_title)
        axTcsc.set_yscale('log')
        axTdcc.set_yscale('log')
        if saveQ:
            figTcsc.savefig(plot_drty+"pop"+plot_save_suffix+".png",dpi=300)
                #replace("pop"+plot_save_suffix+".png",plot_drty+"pop"+plot_save_suffix+".png")

        figTs, axTs = plt.subplots(figsize=(8,8));
        for nn,eln in enumerate(Nsched_vec):
            # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
            # s_lgd_1 = "Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
            for tt, elt in enumerate(xi1_vec):
                xi1 = elt;
                s_xi1 = '\\'+str(elt)+'\\';
                for uu, elu in enumerate(xi2_vec):
                    xi2 = elu;
                    s_xi2 = '\\'+str(elu)+'\\';
                    for ww, elb in enumerate(s_rho_vec):
                        beta = elb;
                        for vv, el in enumerate(sig_vec):
                            s_sig = '(σ='+ str(el) +')';
                            
                            st2_stop = np.flatnonzero(data2[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
                            
                            # if log_yscaleQ:
                            #     axTs.set_yscale('log')
                            axTs.plot(data2[nn][tt][uu][ww][vv][lll_load][0,:st2_stop],data2[nn][tt][uu][ww][vv][lll_load][3,:st2_stop],'--',color=COL2.get_rgb(nn),label='mu frac'+s_lgd_0+s_sig+'_'+str(lll_load))
                            if showNonReprogQ:
                                st1_stop = np.flatnonzero(data1[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
                                axTs.plot(data1[nn][tt][uu][ww][vv][lll_load][0,:st1_stop]+EOT_diff,data1[nn][tt][uu][ww][vv][lll_load][3,:st1_stop],'--',color=COL1.get_rgb(nn),label='mu frac'+s_lgd_1+s_sig+'_'+str(lll_load))

        st0_stop = np.flatnonzero(dataN[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
        axTs.plot(dataN[nn][tt][uu][ww][vv][lll_load][0,:st0_stop],dataN[nn][tt][uu][ww][vv][lll_load][3,:st0_stop],'g--',label='mu frac, no treatment_'+str(lll_load))
        axTs.legend(fancybox=True, framealpha=0.5);
                
        axTs.set_ylabel("mu frac"+log_str)
        axTs.set_xlabel("time (days)")
        axTs.set_title(s_title)
            # if saveQ:
            #     figTs.savefig(plot_drty+"pop"+plot_save_suffix+".png",dpi=300)

        figTsr, axTsr = plt.subplots(figsize=(8,8));
        for nn,eln in enumerate(Nsched_vec):
            # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
            # s_lgd_1 = "Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
            for tt, elt in enumerate(xi1_vec):
                xi1 = elt;
                s_xi1 = '\\'+str(elt)+'\\';
                for uu, elu in enumerate(xi2_vec):
                    xi2 = elu;
                    s_xi2 = '\\'+str(elu)+'\\';
                    for ww, elb in enumerate(rho_vec):
                        beta = elb;
                        for vv, el in enumerate(sig_vec):
                            s_sig = '(σ='+ str(el) +')';
                            
                            st2_stop = np.flatnonzero(data2[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
                            axTsr.plot(data2[nn][tt][uu][ww][vv][lll_load][0,:st2_stop],beta * data2[nn][tt][uu][ww][vv][lll_load][3,:st2_stop] * data2[nn][tt][uu][ww][vv][lll_load][2,:st2_stop] * total_cell_num,'--',color=COL2.get_rgb(nn),label='mu rate'+s_lgd_0+s_sig+'_'+str(lll_load))
                            axTsr.set_ylabel("mu rate"+log_str)
                            axTsr.set_xlabel("time (days)")
                            axTsr.set_title(s_title)
                            if showNonReprogQ:
                                st1_stop = np.flatnonzero(data1[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
                                axTsr.plot(data1[nn][tt][uu][ww][vv][lll_load][0,:st1_stop]+EOT_diff,beta * data1[nn][tt][uu][ww][vv][lll_load][3,:st1_stop] * data1[nn][tt][uu][ww][vv][lll_load][2,:st1_stop] * total_cell_num,'--',color=COL1.get_rgb(nn),label='mu rate'+s_lgd_1+s_sig+'_'+str(lll_load))
        
        st0_stop = np.flatnonzero(dataN[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
        axTsr.plot(dataN[nn][tt][uu][ww][vv][lll_load][0,:st0_stop],beta * dataN[nn][tt][uu][ww][vv][lll_load][3,:st0_stop] * dataN[nn][tt][uu][ww][vv][lll_load][2,:st0_stop] * total_cell_num,'g--',label='mu rate, no treatment_'+str(lll_load))
        axTsr.legend(fancybox=True, framealpha=0.5);
        # if saveQ:
        #     figTsr.savefig(plot_drty+"pop"+plot_save_suffix+".png",dpi=300)
    
    

        figTo, axTo = plt.subplots(figsize=(8,8));
        for nn,eln in enumerate(Nsched_vec):
            # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
            # s_lgd_1 = "Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
            for tt, elt in enumerate(xi1_vec):
                xi1 = elt;
                s_xi1 = '\\'+str(elt)+'\\';
                for uu, elu in enumerate(xi2_vec):
                    xi2 = elu;
                    s_xi2 = '\\'+str(elu)+'\\';
                    for ww, elb in enumerate(rho_vec):
                        beta = elb;
                        for vv, el in enumerate(sig_vec):
                            s_sig = '(σ='+ str(el) +')';
                            
                            st2_stop = np.flatnonzero(data2[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
                            
                            if log_yscaleQ:
                                axTo.set_yscale('log')
                            axTo.plot(data2[nn][tt][uu][ww][vv][lll_load][0,:st2_stop],(data2[nn][tt][uu][ww][vv][lll_load][1,:st2_stop]+data2[nn][tt][uu][ww][vv][lll_load][2,:st2_stop]) * total_cell_num,'--',color=COL2.get_rgb(nn),label='total '+s_lgd_0+s_sig+'_'+str(lll_load))
                            axTo.set_ylabel("Total Cell Number"+log_str)
                            axTo.set_xlabel("time (days)")
                            axTo.set_title(s_title)
                            if showNonReprogQ:
                                st1_stop = np.flatnonzero(data1[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
                                axTo.plot(data1[nn][tt][uu][ww][vv][lll_load][0,:st1_stop]+EOT_diff,(data1[nn][tt][uu][ww][vv][lll_load][1,:st1_stop]+data1[nn][tt][uu][ww][vv][lll_load][2,:st1_stop]) * total_cell_num,'--',color=COL1.get_rgb(nn),label='total '+s_lgd_1+s_sig+'_'+str(lll_load))

        st0_stop = np.flatnonzero(dataN[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
        axTo.plot(dataN[nn][tt][uu][ww][vv][lll_load][0,:st0_stop],(dataN[nn][tt][uu][ww][vv][lll_load][1,:st0_stop]+dataN[nn][tt][uu][ww][vv][lll_load][2,:st0_stop]) * total_cell_num,color[2]+'--',label='total, no treatment_'+str(lll_load))
        axTo.legend(fancybox=True, framealpha=0.5);
        # if saveQ:
        #     figTo.savefig("total_pop"+plot_save_suffix+".png",dpi=300)
        #     # replace(base_dirGD+"total_pop"+plot_save_suffix+".png",plot_drty+"total_pop"+plot_save_suffix+".png")
            
    # TODO: maybe add the other ones too?

        figC, axC = plt.subplots(figsize=(8,8));
        for nn,eln in enumerate(Nsched_vec):
            # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
            # s_lgd_1 = "Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
            for tt, elt in enumerate(xi1_vec):
                xi1 = elt;
                s_xi1 = '\\'+str(elt)+'\\';
                for uu, elu in enumerate(xi2_vec):
                    xi2 = elu;
                    s_xi2 = '\\'+str(elu)+'\\';
                    for ww, elb in enumerate(rho_vec):
                        beta = elb;
                        for vv, el in enumerate(sig_vec):
                            s_sig = '(σ='+ str(el) +')';
                            
                            csc_p = [];
                            csc_p += [data2[nn][tt][uu][ww][vv][lll_load][1]/(data2[nn][tt][uu][ww][vv][lll_load][1] + data2[nn][tt][uu][ww][vv][lll_load][2])];
                            csc_p += [data1[nn][tt][uu][ww][vv][lll_load][1]/(data1[nn][tt][uu][ww][vv][lll_load][1] + data1[nn][tt][uu][ww][vv][lll_load][2])];
                            
                            st2_stop = np.flatnonzero(data2[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
                            
                            
                            # if log_yscaleQ:
                            #     axC.set_yscale('log')
                            axC.plot(data2[nn][tt][uu][ww][vv][lll_load][0,:st2_stop],csc_p[0][:st2_stop],'-',color=COL2.get_rgb(nn),label='CSC frac'+s_lgd_0+s_sig+'_'+str(lll_load))
            
                            if showNonReprogQ:
                                st1_stop = np.flatnonzero(data1[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
                                axC.plot(data1[nn][tt][uu][ww][vv][lll_load][0,:st1_stop]+EOT_diff,csc_p[1][:st1_stop],'-',color=COL1.get_rgb(nn),label='CSC frac'+s_lgd_1+s_sig+'_'+str(lll_load))
            # replace(base_dirGD+"csc"+plot_save_suffix+".png",plot_drty+"csc"+plot_save_suffix+".png")
        csc_p += [dataN[nn][tt][uu][ww][vv][lll_load][1]/(dataN[nn][tt][uu][ww][vv][lll_load][1] + dataN[nn][tt][uu][ww][vv][lll_load][2])];
        st0_stop = np.flatnonzero(dataN[nn][tt][uu][ww][vv][lll_load][0] < T_stop).max();
        axC.plot(dataN[nn][tt][uu][ww][vv][lll_load][0,:st0_stop],csc_p[2][:st0_stop],color[2]+'-',label='CSC frac, no treatment')
        axC.legend(fancybox=True, framealpha=0.5);
        axC.set_ylabel("CSC frac")
        axC.set_xlabel("time (days)")
        axC.set_title(s_title)
        # if saveQ:
        #     figC.savefig("csc"+plot_save_suffix+".png",dpi=300)
            
else:
    s_sched_title = ";reprog:"
    s_title_tweakable = s_sched_title+s_reprg0;
    
    s_lgd_0 = s_sched0+EOT_str;
    s_lgd_1 = s_sched1;
    COL2 = MplColorHelper('viridis', 0, 11)
    COL1 = MplColorHelper('cool', 0, 11)
    figT, axT = plt.subplots(figsize=(8,8));
    for lll_load in load_temp:
        for vv, el in enumerate(sig_vec):
            s_sig = '('+ str(el) +')';
            s_title = "(l:"+f"{float(l_vec[lll_load]):.3}"+";h1:"+f"{float(h1_vec[lll_load]):.3}"+";h2:"+f"{float(h2_vec[lll_load]):.3}"+";hd:"+f"{hd:.3}"+s_title_tweakable+")";
        
            st2_stop = np.flatnonzero(data2[vv][lll_load][0] < T_stop).max();
            st1_stop = np.flatnonzero(data1[vv][lll_load][0] < T_stop).max();
            # st0_stop = np.flatnonzero(dataN[lll_load][0] < T_stop).max();
            # if log_yscaleQ:
            #     axT.set_yscale('log')
            axT.plot(data2[vv][lll_load][0,:st2_stop],data2[vv][lll_load][1,:st2_stop] * total_cell_num,'--',color=COL2.get_rgb(lll_load),label='CSC '+s_lgd_0+s_sig+'_'+str(lll_load))
            axT.plot(data2[vv][lll_load][0,:st2_stop],data2[vv][lll_load][2,:st2_stop] * total_cell_num,'-',color=COL2.get_rgb(lll_load),label='DCC '+s_lgd_0+s_sig+'_'+str(lll_load))
            axT.plot(data1[vv][lll_load][0,:st1_stop],data1[vv][lll_load][1,:st1_stop] * total_cell_num,color[1]+'--',label='CSC '+s_lgd_1)
            axT.plot(data1[vv][lll_load][0,:st1_stop],data1[vv][lll_load][2,:st1_stop] * total_cell_num,color[1]+'-',label='DCC '+s_lgd_1)
            # axT.plot(dataN[lll_load][0,:st0_stop],dataN[lll_load][1,:st0_stop] * total_cell_num,color[2]+'--',label='CSC, no treatment')
            # axT.plot(dataN[lll_load][0,:st0_stop],dataN[lll_load][2,:st0_stop] * total_cell_num,color[2]+'-',label='DCC, no treatment')
            axT.legend(fancybox=True, framealpha=0.5);
            axT.set_ylabel("Cell Number"+log_str)
            axT.set_xlabel("time (days)")
            axT.set_title(s_title)
        axT.set_yscale('log')
        if saveQ:
            figT.savefig(plot_drty+"pop"+plot_save_suffix+".png",dpi=300)
            #replace("pop"+plot_save_suffix+".png",plot_drty+"pop"+plot_save_suffix+".png")
    figTs, axTs = plt.subplots(figsize=(8,8));
    for lll_load in load_temp:
        for vv, el in enumerate(sig_vec):
            s_sig = '('+ str(el) +')';
            s_title = "(l:"+f"{float(l_vec[lll_load]):.3}"+";h1:"+f"{float(h1_vec[lll_load]):.3}"+";h2:"+f"{float(h2_vec[lll_load]):.3}"+";hd:"+f"{hd:.3}"+s_title_tweakable+")";
    
            st2_stop = np.flatnonzero(data2[vv][lll_load][0] < T_stop).max();
            st1_stop = np.flatnonzero(data1[vv][lll_load][0] < T_stop).max();
            # st0_stop = np.flatnonzero(dataN[lll_load][0] < T_stop).max();
            # if log_yscaleQ:
            #     axT.set_yscale('log')
            axTs.plot(data2[vv][lll_load][0,:st2_stop],data2[vv][lll_load][3,:st2_stop],'--',color=COL2.get_rgb(lll_load),label='mu '+s_lgd_0+s_sig+'_'+str(lll_load))
            # axTs.plot(data1[vv][lll_load][0,:st1_stop],data1[vv][lll_load][1,:st1_stop] * total_cell_num,color[1]+'--',label='CSC '+s_lgd_1)
            # axTs.plot(data1[vv][lll_load][0,:st1_stop],data1[vv][lll_load][2,:st1_stop] * total_cell_num,color[1]+'-',label='DCC '+s_lgd_1)
            # axT.plot(dataN[lll_load][0,:st0_stop],dataN[lll_load][1,:st0_stop] * total_cell_num,color[2]+'--',label='CSC, no treatment')
            # axT.plot(dataN[lll_load][0,:st0_stop],dataN[lll_load][2,:st0_stop] * total_cell_num,color[2]+'-',label='DCC, no treatment')
            axTs.legend(fancybox=True, framealpha=0.5);
            axTs.set_ylabel("mu"+log_str)
            axTs.set_xlabel("time (days)")
            axTs.set_title(s_title)
            if saveQ:
                figTs.savefig(plot_drty+"pop"+plot_save_suffix+".png",dpi=300)
    
    
    figTo, axTo = plt.subplots(figsize=(8,8));
    for lll_load in load_temp:
        for vv, el in enumerate(sig_vec):
            s_sig = '('+ str(el) +')';
            s_title = "(l:"+f"{float(l_vec[lll_load]):.3}"+";h1:"+f"{float(h1_vec[lll_load]):.3}"+";h2:"+f"{float(h2_vec[lll_load]):.3}"+";hd:"+f"{hd:.3}"+s_title_tweakable+")";
    
            st2_stop = np.flatnonzero(data2[vv][lll_load][0] < T_stop).max();
            st1_stop = np.flatnonzero(data1[vv][lll_load][0] < T_stop).max();
            # st0_stop = np.flatnonzero(dataN[lll_load][0] < T_stop).max();
            if log_yscaleQ:
                axTo.set_yscale('log')
            axTo.plot(data2[vv][lll_load][0,:st2_stop],(data2[vv][lll_load][1,:st2_stop]+data2[vv][lll_load][2,:st2_stop]) * total_cell_num,'--',color=COL2.get_rgb(lll_load),label='total '+s_lgd_0+s_sig+'_'+str(lll_load))
            axTo.plot(data1[vv][lll_load][0,:st1_stop],(data1[vv][lll_load][1,:st1_stop]+data1[vv][lll_load][2,:st1_stop]) * total_cell_num,color[1]+'--',label='total '+s_lgd_1)
            #axTo.plot(dataN[lll_load][0,:st0_stop],(dataN[lll_load][1,:st0_stop]+dataN[lll_load][2,:st0_stop]) * total_cell_num,color[2]+'--',label='total, no treatment')
            axTo.legend(fancybox=True, framealpha=0.5);
            axTo.set_ylabel("Total Cell Number"+log_str)
            axTo.set_xlabel("time (days)")
            axTo.set_title(s_title)
            if saveQ:
                figTo.savefig("total_pop"+plot_save_suffix+".png",dpi=300)
                # replace(base_dirGD+"total_pop"+plot_save_suffix+".png",plot_drty+"total_pop"+plot_save_suffix+".png")
    
    # TODO: maybe add the other ones too?
    figC, axC = plt.subplots(figsize=(8,8));
    for lll_load in load_temp:
        for vv, el in enumerate(sig_vec):
            s_sig = '('+ str(el) +')';
            s_title = "(l:"+f"{float(l_vec[lll_load]):.3}"+";h1:"+f"{float(h1_vec[lll_load]):.3}"+";h2:"+f"{float(h2_vec[lll_load]):.3}"+";hd:"+f"{hd:.3}"+s_title_tweakable+")";
        
            csc_p = [];
            csc_p += [data2[vv][lll_load][1]/(data2[vv][lll_load][1] + data2[vv][lll_load][2])];
            csc_p += [data1[vv][lll_load][1]/(data1[vv][lll_load][1] + data1[vv][lll_load][2])];
            #csc_p += [dataN[lll_load][1]/(dataN[lll_load][1] + dataN[lll_load][2])];
            st2_stop = np.flatnonzero(data2[vv][lll_load][0] < T_stop).max();
            st1_stop = np.flatnonzero(data1[vv][lll_load][0] < T_stop).max();
            #st0_stop = np.flatnonzero(dataN[lll_load][0] < T_stop).max();
            
            if log_yscaleQ:
                axC.set_yscale('log')
            axC.plot(data2[vv][lll_load][0,:st2_stop],csc_p[0][:st2_stop],'-',color=COL2.get_rgb(lll_load),label='CSC '+s_lgd_0+s_sig+'_'+str(lll_load))
            axC.plot(data1[vv][lll_load][0,:st1_stop],csc_p[1][:st1_stop],color[1]+'-',label='CSC '+s_lgd_1)
            #axC.plot(dataN[lll_load][0,:st0_stop],csc_p[2][:st0_stop],color[2]+'-',label='CSC, no treatment')
            axC.legend(fancybox=True, framealpha=0.5);
            axC.set_ylabel("CSC frac")
            axC.set_xlabel("time (days)")
            axC.set_title(s_title)
            if saveQ:
                figC.savefig("csc"+plot_save_suffix+".png",dpi=300)

if compBaselineQ:
    s_fr_0 = '_'+str(Fracs[0])+'_Fracs_';
    s_fr_1 = '_'+str(Fracs[1])+'_Fracs_';#'_'+str(Fracs[1])+'_Fracs_';
    base_dir = "G:\\My Drive\\a PhD Projects\\GBM Modeling\\python scripts\\data\\k2_model\\"
    caseC = "baseline_beta\\5_Nov"
    base_dirC = base_dir+"conventional_BED\\"+caseC+"\\death_val_of_kim\\corrected_reprog";
    data2C = np.ones((len(rho_vec),len(sig_vec),12)).tolist();
    data1C = np.ones((len(rho_vec),len(sig_vec),12)).tolist();
    dataNC = np.ones((len(rho_vec),len(sig_vec),12)).tolist();
    T_stop = 500;
    rho_vec = [0];
    sig_vec = [0];
    for lll_load in load_temp:
        for ww, elb in enumerate(rho_vec):
            for vv, el in enumerate(sig_vec):
                plot_str2 = "TU_"+str(C[0])+'_reprog'+str(Doses[0])+'_'+comp_str+s_fr_0+str(T_stop)+'days_'+deathFdbk_str+"_"+str(lll_load);
                plot_str1 = "TU_"+str(C[1])+'_reprog'+str(Doses[1])+'_'+comp_str+s_fr_1+str(T_stop)+'days_'+deathFdbk_str+"_"+str(lll_load);
                data2C[ww][vv][lll_load] = np.loadtxt(base_dirC + "\\0\\0" + "\\" + plot_str2 + ".txt");
                data1C[ww][vv][lll_load] = np.loadtxt(base_dirC + "\\0\\0" + "\\" + plot_str1 + ".txt");
                plotn_str2 = "TU_none_"+str(C[0])+'_reprog'+str(Doses[0])+'_'+comp_str+s_fr_0+str(T_stop)+'days_'+deathFdbk_str+"_"+str(lll_load); 
    figVec = [figTcsc, figTdcc, figTs, figTo, figC, figTsr];
    axVec = [axTcsc, axTdcc, axTs, axTo, axC, axTsr];
    sig_vec = [1];
    for ww, elb in enumerate(rho_vec):
        for vv, el in enumerate(sig_vec):
            stC_stop = np.flatnonzero(data2C[ww][vv][lll_load][0] < T_stop).max();
            s_sig = '('+ str(el) +')';
            s_null = ', RT only no mu'
            axVec[0].plot(data2C[ww][vv][lll_load][0,:stC_stop],data2C[ww][vv][lll_load][1,:stC_stop] * total_cell_num,'c--',label='CSC'+s_null+s_lgd_0+s_sig+'_'+str(lll_load))
            axVec[1].plot(data2C[ww][vv][lll_load][0,:stC_stop],data2C[ww][vv][lll_load][2,:stC_stop] * total_cell_num,'c-',label='DCC'+s_null+s_lgd_0+s_sig+'_'+str(lll_load))
            axVec[2].plot(data2C[ww][vv][lll_load][0,:stC_stop],data2C[ww][vv][lll_load][3,:stC_stop],'c:',label='mu'+s_null+s_sig+'_'+str(lll_load))
            axVec[3].plot(data2C[ww][vv][lll_load][0,:stC_stop],(data2C[ww][vv][lll_load][1,:stC_stop]+data2C[ww][vv][lll_load][2,:stC_stop]) * total_cell_num,'c--',label='total'+s_null+s_lgd_0+s_sig+'_'+str(lll_load))
            csc_p_0 = data2C[ww][vv][lll_load][1]/(data2C[ww][vv][lll_load][1] + data2C[ww][vv][lll_load][2])
            axVec[4].plot(data2C[ww][vv][lll_load][0,:stC_stop],csc_p_0[:stC_stop],'c-',label='CSC'+s_null+s_lgd_0+s_sig+'_'+str(lll_load))
            axVec[5].plot(data2C[ww][vv][lll_load][0,:stC_stop],beta * data2C[ww][vv][lll_load][3,:stC_stop] * data2C[ww][vv][lll_load][2,:stC_stop] * total_cell_num,'c-',label='mu rate'+s_null+s_lgd_0+s_sig+'_'+str(lll_load))
    #[axVec[i].legend() for i in [0,1,2,3,4,5]]
    [axVec[i].set_xlim(0,T_stop+100) for i in [0,1,2,3,4,5]]
# # only relevant for feedback regime 6, model 0 or 1.
# eigs2 = np.array(list(map(lambda a,b: calc_eigs6(a,b), data2[vv][lll_load][1,:],data2[vv][lll_load][2,:])))
# eigs1 = np.array(list(map(lambda a,b: calc_eigs6(a,b), data1[vv][lll_load][1,:],data1[vv][lll_load][2,:])))
# #eigs0 = np.array(list(map(lambda a,b: calc_eigs6(a,b), dataN[lll_load][1,:],dataN[lll_load][2,:])))

# for lll_load in load_temp:
#     figE, axE = plt.subplots(figsize=(8,8));
#     for vv, el in enumerate(sig_vec):
#         s_sig = '('+ str(el[0])+','+str(el[1])+')';
    
#         s_title = "(l:"+f"{float(l_vec[lll_load]):.3}"+";h1:"+f"{float(h1_vec[lll_load]):.3}"+";h2:"+f"{float(h2_vec[lll_load]):.3}"+";hd:"+f"{hd:.3}"+s_title_tweakable+")";
        
        
#         st2_stop = np.flatnonzero(data2[vv][lll_load][0] < T_stop).max();
#         st1_stop = np.flatnonzero(data1[vv][lll_load][0] < T_stop).max();
#         # st0_stop = np.flatnonzero(dataN[lll_load][0] < T_stop).max();
    
#         axE.plot(data2[vv][lll_load][0,:st2_stop:50],eigs2[:st2_stop:50,0],'_',color=COL2.get_rgb(vv),label='eig0 '+s_lgd_0+s_sig+'_'+str(lll_load))
#         axE.plot(data2[vv][lll_load][0,:st2_stop:50],eigs2[:st2_stop:50,1],'|',color=COL2.get_rgb(vv),label='eig1 '+s_lgd_0+s_sig+'_'+str(lll_load))
#         axE.plot(data1[vv][lll_load][0,:st1_stop:50],eigs1[:st1_stop:50,0],color[1]+'_',label='eig0 '+s_lgd_1)
#         axE.plot(data1[vv][lll_load][0,:st1_stop:50],eigs1[:st1_stop:50,1],color[1]+'|',label='eig1 '+s_lgd_1)
#         # axE.plot(dataN[lll_load][0,:st0_stop:2],eigs1[:st0_stop:2,0],color[2]+'_',label='eig0, no treatment')
#         # axE.plot(dataN[lll_load][0,:st0_stop:2],eigs1[:st0_stop:2,1],color[2]+'|',label='eig1, no treatment')
#         axE.legend(fancybox=True, framealpha=0.5);
#         axE.set_ylabel("eigs per sim")
#         axE.set_xlabel("time (days)")
#         #axE.set_ylim(-.02,.001)
#         # axE.set_yscale('log')
#         axE.set_title(s_title)

#for i, c in enumerate(cases):
# data_drty = base_dirGD+"\\data\\"+base_model_name+"\\"+model_suffix+"\\"+case+"\\"+date_data_dir+ deathVal_dir + sub_drty;
#G:\\My Drive
#C:\\Users\\jhvo9\\Google Drive (vojh1@uci.edu)

# use_srvQ = True;
# if use_srvQ:
#     sig_vec = [[0,0],[3.6,0.05],[0.05, 3.6]];
#     s_srv_vec = np.arange(len(sig_vec)).tolist();
#     case = 'mu_test3\\26_Oct';
#     for vv, el in enumerate(sig_vec):
#         s_srv_vec[vv] = '\\' + str(el[0]) + '\\' + str(el[1]) + '\\';
# else:
#     sig_vec=[[3.6,0.05]];
#     s_srv_vec = [''];
#     case = 'full_fdbk_gains\\26_Oct';#'full_fdbk_gains\\13_Oct'