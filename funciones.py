# -*- coding: utf-8 -*-
"""
Created on Fri May  7 17:36:25 2021

@author: jhvo9
"""
import numpy as np
import scipy.integrate as integrate

def dU_dt(U,t, r1, r2, d, p, h1, h2, hd, z, l, n, sig, mu_bar, chi):
    #function dU = stem_ODE_feedback(t, U, r1, r2, d, p, h, hd, z, l, n, sig, mu_bar, chi)
    return np.array([(2*p/(1+l*U[1]**n)-1)*r1/(1+h1*U[1]**n)*U[0] + mu_bar * (chi * U[2] / ( 1 + chi * U[2] )) * U[1], 
                     2*(1-p/(1+l*U[1]**n))*r1/(1+h1*U[1]**n)*U[0] + (r2/(1+h2*U[1]**n) - d*hd*U[1]**n/(1+hd*U[1]**n) - mu_bar * chi * U[2] / ( 1 + chi * U[2] )) * U[1],
                    0
                    ])

def radiotherapy_kim(U, LQ_para, surv_vec):
    def fdbk(control, surv):
        val = 1/(1+control*surv);
        return val
    #UNTITLED2 Summary of this function goes here
    #   Detailed explanation goes here
    u, v, s = U[:,-1];
    [a1, b1, a2, b2, c, D] = LQ_para;
    [cont_p_a, cont_p_b, compt_mult, srvn_csc, srvn_dcc, cont_c] = surv_vec;
    # compt_mult tunes fold-change of stem cell feedback ratio based on diff.
    # cell feedback ratio
    SF_U =  np.exp(-a1*fdbk(cont_p_a*compt_mult, s)*D-b1*fdbk(cont_p_b*compt_mult, s)*D**2);
    SF_V = np.exp(-a2*fdbk(cont_p_a, s)*D-b2*fdbk(cont_p_b, s)*D**2);
    u_new = u*SF_U + c*D*v;
    v_new = (v*SF_V - c*D*v); # apply RT; assumes death occurs on non-reprogrammed cells. is this biologically valid?
    # max(v*exp(-a2*D-b2*D^2) - c*v*D,0)
    s_new = s + srvn_csc * (u-u*SF_U) + srvn_dcc * (v-v_new); 
    # v_new is used instead of SF_U because radiotherapy causes de-dif
    return [u_new,v_new, s_new,SF_U, SF_V]

def radiotherapy(U, LQ_para, surv_vec):
    def fdbk(control, surv):
        val = 1/(1+control*surv);
        return val
    u, v, s = U[:,-1];
    [a1, b1, a2, b2, c, D] = LQ_para;
    [cont_p_a, cont_p_b, compt_mult, srvn_csc, srvn_dcc, cont_c] = surv_vec;
    # compt_mult tunes fold-change of stem cell feedback ratio based on diff.
    # cell feedback ratio
    SF_U =  np.exp(-a1*fdbk(cont_p_a*compt_mult, s)*D-b1*fdbk(cont_p_b*compt_mult, s)*D**2);
    SF_V = np.exp(-a2*fdbk(cont_p_a, s)*D-b2*fdbk(cont_p_b, s)*D**2);
    u_new = u*SF_U + min(1,c*D*fdbk(cont_c, s))*v;
    v_new = max(0,(SF_V - min(1,c*D*fdbk(cont_c, s)))*v); # apply RT; assumes death occurs on non-reprogrammed cells. is this biologically valid?
    # max(v*exp(-a2*D-b2*D^2) - c*v*D,0)
    s_new = s + srvn_csc * (u-u*SF_U) + srvn_dcc * (v-v_new); 
    # v_new is used instead of SF_U because radiotherapy causes de-dif
    return [u_new,v_new, s_new,SF_U, SF_V]



def dynamics(para_values, sim_values):
    model0Q, model1Q, model2Q, kimReprogQ, total_cell_num, treat_days, srv_start, LQ_param, total_start_frac, sc_start, sim_resume_days, surv_vec, time_pts1, time_pts2 = sim_values;
    r1, r2, d, p, h1, h2, hd, z, l, n, sig, mu_bar, chi = para_values;
    tc_start = total_start_frac-sc_start;
    #print("data loaded")
    # Defining treatment days and ODE simulation start time after each fraction  
    # ODE simulation before first fraction of RT
    # With treatment 
    U0 = [sc_start, tc_start, srv_start]; 
    T = np.linspace(0, treat_days[0], time_pts1);
    #print("initial growth evaluated")
    U = integrate.odeint(dU_dt, U0, T, args=(r1, r2, d, p, h1, h2, hd, z, l, n, sig, mu_bar, chi)).T
    u_sc, u_dc, u_srv = U
    if model1Q or model2Q:
        u_srv = np.multiply(u_srv, np.exp(-sig * (T))); # using an integrating factor to bypass stiffness
    
    ## pre-therapy growth dynamics and plotting

    ## radiotherapy dynamics and plotting
    #print("repeating")
    for i in range(len(sim_resume_days)):
        
    ####### with stem cell #########
        if kimReprogQ:
            [u_new,v_new, s_new,SF_U, SF_V] = radiotherapy_kim(U, LQ_param, surv_vec);
        else:
            [u_new,v_new, s_new,SF_U, SF_V] = radiotherapy(U, LQ_param, surv_vec);
        
        #print(i,"reprogramming done")
        # integrate.odeint()
        T_int = np.linspace(sim_resume_days[i], treat_days[i+1],int(np.round(time_pts2*(-sim_resume_days[i] + treat_days[i+1]))));
        U0_int = [u_new, v_new, s_new];
        
        U_new = integrate.odeint(dU_dt, U0_int, T_int, args=(r1, r2, d, p, h1, h2, hd, z, l, n, sig, mu_bar, chi)).T
        if model1Q or model2Q:
            U_new[2,:] = np.multiply(U[2,:],np.exp(-sig * (T_int-treat_days[i]))); # using an integrating factor to bypass stiffness
        u_sc, u_dc, u_srv = U_new;
        U = np.hstack((U, U_new))
        T = np.concatenate((T, T_int))
    #print("done")
    T_none = np.linspace(0, treat_days[-1],len(T))
    U_none = integrate.odeint(dU_dt, U0, T_none, args=(r1, r2, d, p, h1, h2, hd, z, l, n, sig, mu_bar, chi)).T

    return U, T, U_none, T_none

def parameter_setup(switch_vec, misc_pars):
    #weak_feedbackQ
    subSelectQ, srvQ, compDosesQ, deathFdbkQ, c_dep_sQ, kimReprogQ, kimDeathValQ, kimICQ, model0Q, model1Q, model2Q = switch_vec;
    DT, post_therapy_end, pwr, h, l_vec, ss = misc_pars;
    
    # Initial Condition and rate
    if kimICQ:
        total_start_frac = 0.0005/64; # Victoria: 0.2/64; Nayeon: 0.0005/64
    else:
        total_start_frac = 0.2/64;
    
    # death rate of DCC
    if kimDeathValQ:
        d  = 1/5 * np.log(2)/DT; 
        hd_str_mod = '_kim_val'
    else:
        d  = np.log(2)/DT; 
        hd_str_mod = '_yu_val';
        
    # compare doses or compare reprogramming rates?
    if compDosesQ:
        v_suffix = ["2x30","1x67"]; # 4 x 13
        v_suffix_save_A = "dose_comp";
        v_suffix_save = "dose_comp_"+str(post_therapy_end);#"dose_comp_0_reprog";
        Doses = [2, 1]; # dose fraction sizes in Gy
        Frac = [30, 67]; 
        C = [5.196*10 ** (-pwr)];
        fg_subtitle = ['Dediff Rate:',str(C[0])];
    else:
        v_suffix = ["w/o reprog","w/ reprog"]; 
        v_suffix_save_A = "dediff_comp";
        v_suffix_save = "dediff_comp_"+str(post_therapy_end);
        Doses = [2]; # dose fraction sizes in Gy
        Frac = [30]; 
        C = [0, 5.196*10**(-pwr)];
        fg_subtitle = 'Dose:'+str(Doses[0])+' Gy, Days:'+str(Frac[0]);
    
    # legacy option of reprogramming or the corrected one?
    if kimReprogQ:
        fg_subtitle = [fg_subtitle, ';original reprog.'];
        reprog_suffix = 'orig_reprog';
    else:
        fg_subtitle = [fg_subtitle, ';improved reprog.'];
        reprog_suffix = 'improved_reprog';
    
    # which kind of death feedback to use? some negative or none?
    if deathFdbkQ:
        hd = 32520.32520325203;#h;
        hd_str = '_death_fdbk';
        hd_suffix = ", and d";
    else:
        hd = 0;
        hd_str = "_no_death_fdbk";
        hd_suffix = '';
    
    # # which kind of feedback onto self-renewal probability? weak or strong?
    # if weak_feedback_Q:
    #     l = 10**(-7); # weak feedback
    #     fdbk_type = 'Weak';
    # else:
    #     l = 10**(3); # strong feedback 
    #     fdbk_type = 'Strong';
    
    # which model type is being used?
    cont_p_a = 0; chi = 0; mu_bar = 0; group_name = 'kim';
    if model0Q == True and model1Q == False and model2Q == False:
        group_name = 'm0';
    elif model0Q == False and model1Q == True and model2Q == False:
        cont_p_a = 10**4;
        group_name = 'm1';
    elif model0Q == False and model1Q == False and model2Q == True:
        cont_p_a = 10**4;
        chi = 10**4;
        mu_bar = C[-1];
        group_name = 'm2';
    else:
        print("this ain't it chief")
    
    # will survivin be added?
    if srvQ:
        #srvn_csc = 3.6; srvn_dcc = 0.05;
        # model is more sensitive to this than to cont_p
        zeta_mult1 = 1;
        zeta_mult2 = 1; #; 3.6, 0.5; 3.6, 5];
    else:
        #srvn_csc = 0; srvn_dcc = 0;
        zeta_mult1 = 0;
        zeta_mult2 = 0;
    
    # will reprogramming rate c depend on survivin?
    cont_c = 0;
    c_fdbk_mod = 'no_c_fdbk';
    if c_dep_sQ:
        cont_c = 10**5;
        c_fdbk_mod = 'yes_c_fdbk';
    
    # will all major kinds of feedback be used?
    if subSelectQ:
        # the point is to force you to pick a value for ss before you can
        # proceed
        rng = [ss];
    else:
        rng = list(range(len(l_vec)));
    # no feedback, only div feedback, only prob feedback (weak, strong), 
    # both div and prob feedback (weak, strong)
    params = [total_start_frac, d, Doses, Frac, C, cont_p_a, chi, mu_bar, hd, zeta_mult1, zeta_mult2, cont_c, rng];
    stgs = [hd_str_mod, v_suffix, v_suffix_save_A, v_suffix_save, fg_subtitle, reprog_suffix, hd_str, hd_suffix, group_name, c_fdbk_mod];
    return [params, stgs]

def data_saver():
    return