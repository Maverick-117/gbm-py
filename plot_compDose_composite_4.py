# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 22:43:58 2021

PLOT dose curves for varying BED for a fixed time.
and
PLOT dose curves for varying time for a fixed BED.

updated to use mu, rho, and sigma
@author: jhvo9
"""

import matplotlib.pyplot as plt;
import numpy as np;
from os import makedirs
from os.path import exists
import matplotlib;
from funciones import calc_EOT, get_EOT_index_pre, calc_BED_Frac
import matplotlib as mpl
from matplotlib import cm

#### SETUP
class MplColorHelper:

  def __init__(self, cmap_name, start_val, stop_val):
    self.cmap_name = cmap_name
    self.cmap = plt.get_cmap(cmap_name)
    self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
    self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

  def get_rgb(self, val):
    return self.scalarMap.to_rgba(val)


font = {'weight' : 'bold',
        'size'   : 13}

matplotlib.rc('font', **font)

total_cell_num = 4/3*np.pi*(1 ** 3)*1e9;
goldenLaptopQ = False;
compDosesQ = True; #true: yes you're comparing doses; false: no, you're comparing reprogramming versus non-reprogramming
saveQ = False;
log_yscaleQ = False;
reprogQ = True;
if reprogQ:
    rp_str = 'reprog';
else:
    rp_str = 'no reprog';
r2 = np.log(2)/3.9;
n = 1; treat_start = 100;

subSelectQ = True; kimReprogQ = False; kimDeathValQ = True; 
l_w = 10**(-7); # weak feedback on prob
l_s = 10**3; # strong feedback on prob
h1 = 10**5; # feedback on css div 
h2 = 10**5; # feedback on dcc div
l_vec =  [0,  0, l_w, l_s, l_w , l_s , 0 ,0 ,l_w,l_w,l_s,l_s];
h1_vec = [0, h1, 0  , 0  , h1  , h1  , h1,0 ,h1 ,0  , h1, 0 ];
h2_vec = [0, h2, 0  , 0  , h2  , h2  , 0 ,h2,0  ,h2 , 0 , h2];
#pre_lll_load = 9;
if subSelectQ:
    #lll_vec = [pre_lll_load];
    #lll_vec = [4,8,9];
    lll_vec = [6];
    #[4,8,9,5,10,11];#[0,1,3,4,8,9,5,10,11];#[4,8,9,5,10,11];
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
comp_conventional_60_30Q = True;

### SETTINGS
base_model_name = 'k2_model'
#case = "34";
model_suffix = "conventional_BED"#
date_data_dir = '1_Dec\\'
date_plot_dir = '1_Dec\\'; 
# 11_Dec, T_stop = 1200 and add some post-EOT times to see the minor but clear evidence that post-1 year dynamics aren't flat
if kimDeathValQ:
    deathVal_dir = 'death_val_of_kim\\';
else:
    deathVal_dir = 'death_val_of_not_kim\\';

if kimReprogQ:
    pre_sub_drty = "kim_reprog\\";
else:
    pre_sub_drty = "corrected_reprog\\";

comp_str = "Gy";
a,b =  np.array([0.17, 0.02]); 
# Doses = [1, 2, 40/15,34/10,5];
# Fracs = [60,30,15   ,10   ,5];
DosesEmp = [1, 2, 40/15,34/10,5];
FracsEmp = [60,30,15   ,10   ,5];
TotalDoseEmp = [int(DosesEmp[i]*FracsEmp[i]) for i in range(len(DosesEmp))];
BEDEmp = [TotalDoseEmp[i] * (1+DosesEmp[i]/8.5) for i in range(len(DosesEmp))];
base_idxs = list(range(1,len(DosesEmp)));
sig_list = sorted(list(map(lambda p:  10 ** (p), np.arange(-2,3).tolist()))+[0]);
rho_list = [0.2];#sorted(list(map(lambda p: 2 ** (p), np.arange(0,-3,-1).tolist())) + [0,2]);
# rename rho to rho
xi1_list = [0];#sorted(list(map(lambda p:  10 ** (p), np.arange(-2,3).tolist()))+[0]);
xi2_list = sorted(list(map(lambda p:  10 ** (p), np.arange(-2,3).tolist()))+[0]);#np.arange(0.01,0.1,0.01);#[.1];#sorted(list(map(lambda p:  10 ** (p), np.arange(0,3).tolist()))+[0])
days_past_treatment = [30,90,180,360];
pick_idx = [1];
T_stop = 600;#int(t_vec.max());

Doses = sorted(np.arange(1,21,dtype=float).tolist() + DosesEmp[2:5]);

#collect time points in the space of parameters we're working with
#results_mtx = np.ones((4,len(Doses),len(days_past_treatment),4)); # BED x Doses x DPT x cell-types

tp = np.ones((len(xi1_list),len(xi2_list),len(rho_list),len(sig_list),len(Doses)));
results_mtx = np.ones((1,len(xi1_list),len(xi2_list),len(rho_list),len(sig_list),len(Doses),len(days_past_treatment),6)); # BED x rho x sigma x Doses x DPT x (cell-types and mu)
for p,ref_idx in enumerate(pick_idx): # per BED (p)

    Fracs = np.array(list(map(lambda d:calc_BED_Frac(a,b,d,BEDEmp[ref_idx]),Doses)));
    print("plot loop for total dose",p,ref_idx,TotalDoseEmp[ref_idx])
    #Doses = sorted(np.arange(1,21,dtype=int).tolist() + [DosesEmp[ref_idx]]); # dose fraction sizes in Gy
    
    if not comp_conventional_60_30Q:
        pulley = np.where(np.array(Doses) == DosesEmp[ref_idx])[0][0];
    else:
        pulley = 1;#np.where(np.array(Doses) == DosesEmp[1])[0][0]
    print(pulley)
    if ref_idx in [2,3]:
        # case = str(TotalDoseEmp[ref_idx]);
        case_str = str(np.round(DosesEmp[ref_idx],decimals=2))+" Gy";
    elif ref_idx in [1]:
        # case = 'reversionAttempt5'
        case_str = "conventional";
    elif ref_idx in [0,4]:
        # case = str(TotalDoseEmp[ref_idx]);
        case_str = str(DosesEmp[ref_idx])+" Gy";
    case = 'debug2'
    # case_str = "conventional";
    C = [ 5.196*10**(-3), 0.0];
    fix_str = str(C[0])+"_reprog";
    comp_dir = "Gy_vs_conv"
    comp_list = Doses;
    unit = ' Gy';
    csr = Doses[::-1];
    color = ['r','k','g','b','m'];
    
        
    if deathFdbkQ:
        hd = 32520.32520325203;#1e5; 
        deathFdbk_str = "_w_death_fdbk";
    else:
        hd = 0.0; 
        deathFdbk_str= "_w_no_death_fdbk";
    ### PLOTTING
    for lll_load in lll_vec:
        total_vec = tp.tolist(); csc_frac_vec = tp.tolist();
        totaln_vec = tp.tolist(); cscn_frac_vec = tp.tolist();
        EOT_vec = tp.tolist();
        time_vec = tp.tolist(); timen_vec = tp.tolist();
        mu_vec = tp.tolist(); mun_vec = tp.tolist();
        # plotting prep per feedback type
        print("evaluating and graphic index:", lll_load)
        l = float(l_vec[lll_load]);
        h1 = float(h1_vec[lll_load]);
        h2 = float(h2_vec[lll_load]);
        probFdbkQ = l > 0;
        divR1FdbkQ = h1 > 0;
        divR2FdbkQ = h2 > 0;
        #directory hunting
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
        for tt, elt in enumerate(xi1_list):
            xi1 = elt;
            s_xi1 = '\\'+str(elt);
            for uu, elu in enumerate(xi2_list):
                xi2 = elu;
                s_xi2 = '\\'+str(elu);
                for ww, elr in enumerate(rho_list):
                    rho = elr;
                    s_rho = '\\'+str(elr);
                    for vv, el in enumerate(sig_list):
                        sig = el;
                        s_sig = '\\'+ str(el) + '\\';
                        s_mu_suffix = s_xi1 + s_xi2 + s_rho + s_sig+ 'Schedule0\\';
                        sub_drty = pre_sub_drty + s_mu_suffix ;
                        data_drty = base_dirGD+"\\a PhD Projects\\GBM Modeling\\python scripts\\data\\"+base_model_name+"\\"+model_suffix+"\\"+case+"\\"+date_data_dir+ deathVal_dir + sub_drty
                        plot_drty = base_dirGD+"\\a PhD Projects\\GBM Modeling\\python scripts\\plots\\"+base_model_name+"\\"+model_suffix+"\\"+case+"\\"+date_plot_dir+deathVal_dir + sub_drty + fdbk_dirs +'\\'+comp_dir+ log_dir
                        if not exists(plot_drty):
                            makedirs(plot_drty);
                        if int(comp_list[pulley]) == comp_list[pulley]:
                            dsg = float(comp_list[pulley]);
                        else:
                            dsg = float(comp_list[pulley]);
                        s_fr_0 = '_'+str(float(Fracs[1]))+'_Fracs_'+str(T_stop)+'days_';
                        s_fr_1 = '';#'_'+str(Fracs[1])+'_Fracs_';
                        plot_str = "TU_"+str(C[0])+"_reprog"+str(dsg)+'_'+comp_str+s_fr_0+deathFdbk_str+"_"+str(lll_load);
                        plotn_str = "TU_"+str(C[1])+"_reprog"+str(dsg)+'_'+comp_str+s_fr_0+deathFdbk_str+"_"+str(lll_load);
                        TU_2Gy = np.loadtxt(data_drty+plot_str+".txt");
                        TU_2Gy_none = np.loadtxt(data_drty+plotn_str+".txt");
                        t2_vec = TU_2Gy[0,:];
                        u2_sc = TU_2Gy[1,:];
                        u2_dc = TU_2Gy[2,:];
                        u2_mu = TU_2Gy[3,:]; # save this directly so that the np.cumsum function can be applied to it later
                        #un2_mu_cumsum = np.cumsum(u2_mu);
                    
                        tn2_vec = TU_2Gy_none[0,:];
                        un2_sc = TU_2Gy_none[1,:];
                        un2_dc = TU_2Gy_none[2,:];
                        un2_mu = TU_2Gy_none[3,:];
                        #un2_mu_cumsum = np.cumsum(un2_mu);
                        total2 = u2_sc + u2_dc;
                        csc_frac2 = u2_sc / total2;
                        totaln2 = un2_sc + un2_dc
                        cscn_frac2 = un2_sc / totaln2;
                    
            
                        for ddd in range(len(comp_list)):#collect per Dose (q)
                            dose_temp = float(comp_list[ddd]);
                            # if int(dose_temp) == dose_temp:
                            #     dose_temp = int(dose_temp);
                            s_fr_0 = '_'+str(Fracs[ddd])+'_Fracs_'+str(T_stop)+'days_';
                            s_fr_1 = '';#'_'+str(Fracs[1])+'_Fracs_
                            plot_str = "TU_"+str(C[0])+"_reprog"+str(dose_temp)+'_'+comp_str+s_fr_0+deathFdbk_str+"_"+str(lll_load);
                            plotn_str = "TU_"+str(C[1])+"_reprog"+str(dose_temp)+'_'+comp_str+s_fr_0+deathFdbk_str+"_"+str(lll_load);
                            TU = np.loadtxt(data_drty+plot_str+".txt");
                            TU_none = np.loadtxt(data_drty+plotn_str+".txt");
                            
                            t_vec = TU[0,:];
                            u_sc = TU[1,:];
                            u_dc = TU[2,:];
                            u_mu = TU[3,:];
                            
                            tn_vec = TU_none[0,:];
                            un_sc = TU_none[1,:];
                            un_dc = TU_none[2,:];
                            un_mu = TU_none[3,:];
                            
                            total_vec[tt][uu][ww][vv][ddd] = u_sc + u_dc;
                            csc_frac_vec[tt][uu][ww][vv][ddd] = u_sc / total_vec[tt][uu][ww][vv][ddd];
                            totaln_vec[tt][uu][ww][vv][ddd] = un_sc + un_dc
                            cscn_frac_vec[tt][uu][ww][vv][ddd] = un_sc / totaln_vec[tt][uu][ww][vv][ddd];
                            time_vec[tt][uu][ww][vv][ddd] = t_vec; 
                            timen_vec[tt][uu][ww][vv][ddd] = tn_vec;
                            # calculate end of treatment day
                            EOT_vec[tt][uu][ww][vv][ddd] = calc_EOT(Fracs[ddd], treat_start);
                            mu_vec[tt][uu][ww][vv][ddd] = u_mu; mun_vec[tt][uu][ww][vv][ddd] = un_mu;
                            # mu_rate_vec[tt][uu][ww][vv][ddd] = 
            
            
                        post_EOT_trackers = np.array([
                            [np.nonzero(np.greater_equal(time_vec[tt][uu][ww][vv][ddd], EOT_vec[tt][uu][ww][vv][ddd]+p))[0][0] for ddd in range(len(Doses))]
                            for p in days_past_treatment]).T;
                        EOT_tracker = [np.nonzero(
                            np.less_equal(time_vec[tt][uu][ww][vv][ddd], EOT_vec[tt][uu][ww][vv][ddd]))[0][-1] for ddd in range(len(Doses))];
                        
                        
                        # total_mtx = np.ones((len(Doses),len(days_past_treatment)));
                        # csc_frac_mtx = np.ones((len(Doses),len(days_past_treatment)));
                        # totaln_mtx = np.ones((len(Doses),len(days_past_treatment)));
                        # csc_fracn_mtx = np.ones((len(Doses),len(days_past_treatment)));
                        # collect cell number counts while either relativizing against (60,30) or against the spawn of the BED
                        if not comp_conventional_60_30Q:
                            # against the spawn of the BED
                            rel_measure_total = total_vec[pulley][post_EOT_trackers[:][:]];
                            rel_measure_totaln = totaln_vec[pulley][post_EOT_trackers[:][:]];
                            rel_measure_csc = csc_frac_vec[pulley][post_EOT_trackers[:][:]];
                            rel_measure_cscn = cscn_frac_vec[pulley][post_EOT_trackers[:][:]];
                            rel_measure_mu = mu_vec[pulley][post_EOT_trackers[:][:]];
                            rel_measure_mun = mun_vec[pulley][post_EOT_trackers[:][:]];
                        else:
                            # against (60,30)
                            data_drty_tmp = base_dirGD+"\\a PhD Projects\\GBM Modeling\\python scripts\\data\\"+base_model_name+"\\"+model_suffix+'\\'+case+'\\'+date_data_dir+ deathVal_dir + pre_sub_drty+ s_mu_suffix
                            s_fr_0 = '_'+str(calc_BED_Frac(a,b,2,BEDEmp[ref_idx]))+'_Fracs_'+str(T_stop)+'days_';        
                            plot_str_tmp = "TU_"+str(C[0])+"_reprog2.0_Gy"+s_fr_0+deathFdbk_str+"_"+str(lll_load);
                            plotn_str_tmp = "TU_"+str(C[1])+"_reprog2.0_Gy"+s_fr_0+deathFdbk_str+"_"+str(lll_load);
                            TU_2Gy = np.loadtxt(data_drty_tmp+plot_str_tmp+".txt");
                            TU_2Gy_none = np.loadtxt(data_drty_tmp+plotn_str_tmp+".txt");
                            t2_vec = TU_2Gy[0,:];
                            u2_sc = TU_2Gy[1,:];
                            u2_dc = TU_2Gy[2,:];
                            u2_mu = TU_2Gy[3,:];
                        
                            tn2_vec = TU_2Gy_none[0,:];
                            un2_sc = TU_2Gy_none[1,:];
                            un2_dc = TU_2Gy_none[2,:];
                            un2_mu = TU_2Gy_none[3,:];
                            total2 = u2_sc + u2_dc;
                            csc_frac2 = u2_sc / total2;
                            totaln2 = un2_sc + un2_dc
                            cscn_frac2 = un2_sc / totaln2;
                            conventional_EOT = calc_EOT(30, treat_start);
                            post_EOT_trackers2 = np.array([np.nonzero(np.greater_equal(t2_vec, conventional_EOT+p))[0][0] for p in days_past_treatment]).T;
                            #post_EOT_trackers2n = np.array([np.nonzero(np.greater_equal(tn2_vec, conventional_EOT+p))[0][0] for p in days_past_treatment]).T;
                            rel_measure_total = np.tile(total2[post_EOT_trackers2[:]],(len(Doses),1));
                            rel_measure_totaln = np.tile(totaln2[post_EOT_trackers2[:]],(len(Doses),1));
                            rel_measure_csc = np.tile(csc_frac2[post_EOT_trackers2[:]],(len(Doses),1));
                            rel_measure_cscn = np.tile(cscn_frac2[post_EOT_trackers2[:]],(len(Doses),1));  
                            rel_measure_mu = np.tile(u2_mu[post_EOT_trackers2[:]],(len(Doses),1));
                            rel_measure_mun = np.tile(un2_mu[post_EOT_trackers2[:]],(len(Doses),1));     
                        #relativizing
                        for ddd in range(len(comp_list)):
                            results_mtx[p,tt,uu,ww,vv,ddd,:,0] = total_vec[tt][uu][ww][vv][ddd][post_EOT_trackers[ddd]]/rel_measure_total[ddd];
                            results_mtx[p,tt,uu,ww,vv,ddd,:,1] = totaln_vec[tt][uu][ww][vv][ddd][post_EOT_trackers[ddd]]/rel_measure_totaln[ddd];
                            results_mtx[p,tt,uu,ww,vv,ddd,:,2] = csc_frac_vec[tt][uu][ww][vv][ddd][post_EOT_trackers[ddd]]/rel_measure_csc[ddd];
                            results_mtx[p,tt,uu,ww,vv,ddd,:,3] = cscn_frac_vec[tt][uu][ww][vv][ddd][post_EOT_trackers[ddd]]/rel_measure_cscn[ddd];
                            # if rel_measure_mu[ddd].max() > 0:
                            #     results_mtx[p,ww,vv,ddd,:,4] = np.cumsum(mu_vec[ww][vv][ddd][post_EOT_trackers[ddd,vv,ww]])/rel_measure_mu[ddd];
                            #     results_mtx[p,ww,vv,ddd,:,5] = np.cumsum(mun_vec[ww][vv][ddd][post_EOT_trackers[ddd,vv,ww]])/rel_measure_mun[ddd];
                            # else:
                            results_mtx[p,tt,uu,ww,vv,ddd,:,4] = np.cumsum(mu_vec[tt][uu][ww][vv][ddd])[post_EOT_trackers[ddd]]
                            results_mtx[p,tt,uu,ww,vv,ddd,:,5] = np.cumsum(mun_vec[tt][uu][ww][vv][ddd])[post_EOT_trackers[ddd]]
            

s_title = "(l:"+f"{l:.3}"+";h1:"+f"{h1:.3}"+";h2:"+f"{h2:.3}"+";hd:"+f"{hd:.3}"+")"
plot_save_suffix = '_'+fix_str+'_'+str(T_stop)+'days_'+deathFdbk_str+'_'+str(lll_load);
dosesLabel = [int(d) if int(d) == d else round(d,2) for d in np.unique(Doses)];
if not comp_conventional_60_30Q:
    tl_suffix = "";
    yl_suffix = "("+case_str+")";
else:
    tl_suffix = "("+case_str+")";
    yl_suffix = "(conventional)";
if reprogQ:
    c1 = 0; c2 = 2;
else:
    c1 = 1; c2 = 3;
#plotting the outcome
# for p,d in enumerate(days_past_treatment):
#     # fixing DPT and letting sigma vary in the graph
#     for tt, elt in enumerate(xi1_list):
#         xi1 = elt;
#         for uu, elu in enumerate(xi2_list):
#             xi2 = elu;
    
#             figT, axT = plt.subplots(figsize=(9,4)); #plotting total size
#             figC, axC = plt.subplots(figsize=(9,4)); #plotting CSC frac

            
#             #for p,d in enumerate(days_past_treatment):
#             COL = MplColorHelper('viridis', 0, len(sig_list)-1)
#             for ww, elr in enumerate(rho_list):
#                 for vv, el in enumerate(sig_list):
#                     for q,b in enumerate(pick_idx): 
#                         # results_mtx: BED x Doses x DPT x cell-types
#                         axT.plot(Doses,results_mtx[q,tt,uu,ww,vv,:,p,c1],'o',color=COL.get_rgb(vv),label=str(el))
#                         axC.plot(Doses,results_mtx[q,tt,uu,ww,vv,:,p,c2],'o',color=COL.get_rgb(vv),label=str(el))
                    
#                         if log_yscaleQ:
#                             axT.set_yscale('log');
#                             axC.set_yscale('log');
#                         axT.set_title("total size; "+rp_str+'; '+str(d)+" days since EOT")
#                         axT.set_ylabel("Total Size(Dose)/Total Size"+yl_suffix,fontsize = 11)
#                         axT.set_xlabel("Dose per Frac (Gy)",fontsize = 11)
#                         legT = axT.legend(fancybox=True, framealpha=0.5, title = "sigma",loc='center left', bbox_to_anchor=(1, 0.5))
#                         axT.set_xticks(dosesLabel)
#                         axT.set_xticklabels(dosesLabel,  rotation = 45, ha="right")
#                         # axT.set(ylim=(0.6,1.6),autoscale_on=False)
#                         legT.set_in_layout(True)
#                         figT.tight_layout()
#                         legC = axC.legend(fancybox=True, framealpha=0.5, title = "sigma",loc='center left', bbox_to_anchor=(1, 0.5))
#                         axC.set_title("csc frac; "+rp_str+'; '+str(d)+" days since EOT")
#                         axC.set_ylabel("CSC Frac(Dose)/CSC Frac"+yl_suffix,fontsize = 11)
#                         axC.set_xlabel("Dose per Frac (Gy)",fontsize = 11)
#                         axC.set_xticks(dosesLabel)
#                         axC.set_xticklabels(dosesLabel, rotation = 45, ha="right")
                        
#                         # axC.set(ylim=(0.8,2.0),autoscale_on=False)
#                         legC.set_in_layout(True)
#                         figC.tight_layout()
#                         if saveQ:
#                             figT.savefig(plot_drty+"\\total_size_dosage_"+str(days_past_treatment[p])+" days"+plot_save_suffix+".png",dpi=300)
#                             figC.savefig(plot_drty+"\\csc_frac_dosage_"+str(days_past_treatment[p])+" days"+plot_save_suffix+".png",dpi=300)

for q,b in enumerate(pick_idx):
    # fixing sigma and letting time vary in the graph
    if b == 1:
        tl_suffix_a = 'conventional)'
    else:
        tl_suffix_a = str(round(DosesEmp[b],2))+' Gy)';
    for tt, elt in enumerate(xi1_list):
        xi1 = elt;
        for uu, elu in enumerate(xi2_list):
            xi2 = elu;
            for ww, elr in enumerate(rho_list):
                for vv, el in enumerate(sig_list):
                    tl_suffix = tl_suffix_a  + 'ξ1:' + str(elt) + ',ξ2:' + str(elu) + ';ρ:'+ str(elr)+ ';σ:' + str(el)
                    figTt, axTt = plt.subplots(figsize=(9,4))
                    figCt, axCt = plt.subplots(figsize=(9,4))
                    COL = MplColorHelper('cool', 0, len(days_past_treatment)-1)
        
                    for p,d in enumerate(days_past_treatment):
                        axTt.plot(Doses,results_mtx[q,tt,uu,ww,vv,:,p,c1],'o',color=COL.get_rgb(p),label=str(d))
                        axCt.plot(Doses,results_mtx[q,tt,uu,ww,vv,:,p,c2],'o',color=COL.get_rgb(p),label=str(d))
                    
                    if log_yscaleQ:
                        # axTt.set_yscale('log');
                        axCt.set_yscale('log');
                    axTt.set_title("total size; "+rp_str+'('+tl_suffix)
                    axTt.set_ylabel("Total Size(Dose)/Total Size"+yl_suffix,fontsize = 11)
                    axTt.set_xlabel("Dose per Frac (Gy)",fontsize = 11)
                    legTt = axTt.legend(fancybox=True, framealpha=0.5, title = "days after EOT",loc='center left', bbox_to_anchor=(1, 0.5))
                    axTt.set_xticks(dosesLabel)
                    axTt.set_xticklabels(dosesLabel,  rotation = 45, ha="right")
                    # axTt.set(ylim=(0.6,1.6),autoscale_on=False)
                    legTt.set_in_layout(True)
                    figTt.tight_layout()
                    legCt = axCt.legend(fancybox=True, framealpha=0.5, title = "days after EOT",loc='center left', bbox_to_anchor=(1, 0.5))
                    axCt.set_title("csc frac; "+rp_str+'('+tl_suffix)
                    axCt.set_ylabel("CSC Frac(Dose)/CSC Frac"+yl_suffix,fontsize = 11)
                    axCt.set_xlabel("Dose per Frac (Gy)",fontsize = 11)
                    axCt.set_xticks(dosesLabel)
                    axCt.set_xticklabels(dosesLabel, rotation = 45, ha="right")
                    # axCt.set(ylim=(0.8,2.0),autoscale_on=False)
                    legCt.set_in_layout(True)
                    figCt.tight_layout()
                    if saveQ:
                        # figTt.savefig(plot_drty+"\\total_size_dosage_("+str(DosesEmp[b])+","+str(FracsEmp[b])+")_"+plot_save_suffix+".png",dpi=300)
                        figCt.savefig(plot_drty+  "\\csc_frac_dosage_("+str(DosesEmp[b])+","+str(FracsEmp[b])+")_"+plot_save_suffix+".png",dpi=300)
    
for q,b in enumerate(pick_idx):
    # fixing sigma and letting time vary in the graph
    if q == 0:
        tl_suffix_a = 'conventional)'
    else:
        tl_suffix_a = str(round(DosesEmp[b],2))+' Gy)';
    for tt, elt in enumerate(xi1_list):
        xi1 = elt;
        for uu, elu in enumerate(xi2_list):
            xi2 = elu;
            for ww, elr in enumerate(rho_list):
                for vv, el in enumerate(sig_list):
                    tl_suffix = tl_suffix_a + 'ξ1:' + str(elt) + ',ξ2:' + str(elu) + ';ρ:'+ str(elr)+ ';σ:' + str(el)
                    figM, axM = plt.subplots(figsize=(11,4))
                    #figMn, axMn = plt.subplots(figsize=(11,4))
                    COL = MplColorHelper('cool', 0, len(days_past_treatment)-1)
        
                    for p,d in enumerate(days_past_treatment):
                        axM.plot(Doses,results_mtx[q,tt,uu,ww,vv,:,p,4],'o',color=COL.get_rgb(p),label=str(d))
                        #axMn.plot(Doses,results_mtx[q,tt,uu,ww,vv,:,p,5],'o',color=COL.get_rgb(p),label=str(d))
                    
                    if log_yscaleQ:
                        axM.set_yscale('log');
                        #axMn.set_yscale('log');
                    axM.set_title("Cum. mu(RT reprog); "+'('+tl_suffix)
                    axM.set_ylabel("Cum. mu (Dose)",fontsize = 11)
                    axM.set_xlabel("Dose per Frac (Gy)",fontsize = 11)
                    legM = axM.legend(fancybox=True, framealpha=0.5, title = "days after EOT",loc='center left', bbox_to_anchor=(1, 0.5))
                    axM.set_xticks(dosesLabel)
                    axM.set_xticklabels(dosesLabel,  rotation = 45, ha="right")
                    # axTt.set(ylim=(0.6,1.6),autoscale_on=False)
                    legM.set_in_layout(True)
                    figM.tight_layout()
                    # legMn = axMn.legend(fancybox=True, framealpha=0.5, title = "days after EOT",loc='center left', bbox_to_anchor=(1, 0.5))
                    # axMn.set_title("Cum. mu(no RT reprog); "+'('+tl_suffix)
                    # axMn.set_ylabel("Cum. mu (Dose)",fontsize = 11)
                    # axMn.set_xlabel("Dose per Frac (Gy)",fontsize = 11)
                    # axMn.set_xticks(dosesLabel)
                    # axMn.set_xticklabels(dosesLabel, rotation = 45, ha="right")
                    # axCt.set(ylim=(0.8,2.0),autoscale_on=False)
                    # legMn.set_in_layout(True)
                    # figMn.tight_layout()
                    if saveQ:
                        figM.savefig(plot_drty+"\\mu_RT_("+str(DosesEmp[b])+","+str(FracsEmp[b])+")_"+plot_save_suffix+".png",dpi=300)
                        # figMn.savefig(plot_drty+  "\\mu_no_RT_("+str(DosesEmp[b])+","+str(FracsEmp[b])+")_"+plot_save_suffix+".png",dpi=300)
