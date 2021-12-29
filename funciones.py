# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 18:56:28 2021

@author: jhvo9
"""
import numpy as np
import scipy.integrate as integrate

def dU_dt(U,t, r1, r2, d, p, h1, h2, hd, z, l, n, sig, mu_bar, rho, xi1, xi2, ddd):
    #function dU = stem_ODE_feedback(t, U, r1, r2, d, p, h, hd, z, l, n, sig, mu_bar, chi)
    return np.array([(2*p/(1+l*U[1]**n)-1)*r1/(1+h1*U[1]**n)*U[0] + rho * xi1*ddd/(1+xi1*ddd) * U[2] * U[1], #- 0*(d+(1.1*r2-d)*hd*U[1]**n/(1+hd*U[1]**n)) * U[0], 
                      2*(1-p/(1+l*U[1]**n))*r1/(1+h1*U[1]**n)*U[0] + U[1] * (r2/(1+h2*U[1]**n) - d-(1.1*r2-d)*hd*U[1]**n/(1+hd*U[1]**n) - rho * xi1*ddd/(1+xi1*ddd) * U[2]),
                      sig * (1/(1+xi2*ddd)) * (mu_bar - U[2])
                    ])
    # make simulations, argue from equations

def calc_BED_Frac(a,b,Dose,baseline_BED=60*(1+2/8.5)):
    # a = .17, b =.02
    # the baseline BED is set to conventional dosage, by default
    return np.floor(baseline_BED/(Dose*(1+Dose/(a/b))))

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
    #radiotherapy_mkiii
    u, v, mu = U[:,-1];
    [a1, b1, a2, b2, c, D] = LQ_para;
    [cont_p_a, cont_p_b, compt_mult, srvn_csc, srvn_dcc, cont_c] = surv_vec;
    SF_U =  np.exp(-a1*D-b1*D**2);
    SF_V = np.exp(-a2*D-b2*D**2);
    mu_new = mu +c * D+ .01428;  # 2 Nov: added "+ .01428" 
    v_new = max(0,(1 - min(1,mu_new))*SF_V*v);
    u_new = u*SF_U + min(1,mu_new)*SF_V*v;
    return [u_new,v_new, mu_new,SF_U, SF_V]

def radiotherapy_mu0(U, LQ_para, surv_vec):
    #radiotherapy mu = 0 check
    u, v, mu = U[:,-1];
    [a1, b1, a2, b2, c, D] = LQ_para;
    [cont_p_a, cont_p_b, compt_mult, srvn_csc, srvn_dcc, cont_c] = surv_vec;
    SF_U =  np.exp(-a1*D-b1*D**2);
    SF_V = np.exp(-a2*D-b2*D**2);
    mu_new = mu +c * D + .01428;  # 2 Nov: added "+ .01428" 
    v_new = max(0,(1 - min(1,mu_new))*SF_V*v);
    u_new = u*SF_U + min(1,mu_new)*SF_V*v;
    mu_new = 0;
    return [u_new,v_new, mu_new,SF_U, SF_V]

def dynamics(para_values, sim_values):
    # dynamics_mkiii
    model0Q, model1Q, model2Q, kimReprogQ, total_cell_num, treat_days, mu_start, LQ_param, total_start_frac, sc_start, sim_resume_days, surv_vec, time_pts1, time_pts2 = sim_values;
    r1, r2, d, p, h1, h2, hd, z, l, n, sig, mu_bar, beta, xi1, xi2, ddd = para_values;
    tc_start = total_start_frac-sc_start;
    # ODE simulation before first fraction of RT
    # With treatment 
    U0 = [sc_start, tc_start, mu_start]; 
    T = np.linspace(0, treat_days[0], time_pts1);
    #print("initial growth evaluated")
    U = integrate.odeint(dU_dt, U0, T, args=(r1, r2, d, p, h1, h2, hd, z, l, n, sig, mu_bar, beta, xi1, xi2, ddd)).T

    
    ## pre-therapy growth dynamics and plotting

    ## radiotherapy dynamics 
    for i in range(len(sim_resume_days)):
        

        if kimReprogQ:
            [u_new,v_new, mu_new,SF_U, SF_V] = radiotherapy_kim(U, LQ_param, surv_vec);
        else:
            [u_new,v_new, mu_new,SF_U, SF_V] = radiotherapy(U, LQ_param, surv_vec);
        
        #print(i,"reprogramming done")
        # integrate.odeint()
        T_int = np.linspace(sim_resume_days[i], treat_days[i+1],int(np.round(time_pts2*(-sim_resume_days[i] + treat_days[i+1]))));
        U0_int = [u_new, v_new, mu_new];

        U_new = integrate.odeint(dU_dt, np.array(U0_int).reshape((3,)), T_int, args=(r1, r2, d, p, h1, h2, hd, z, l, n, sig, mu_bar, beta, xi1, xi2, ddd)).T

        U = np.hstack((U, U_new))
            
        T = np.concatenate((T, T_int))

    #print("done")
    T_none = np.linspace(0, treat_days[-1],time_pts2)
    U_none = integrate.odeint(dU_dt, U0, T_none, args=(r1, r2, d, p, h1, h2, hd, z, l, n, sig, mu_bar, beta, xi1, xi2, ddd)).T
    
    return U, T, U_none, T_none

def calc_EOT(frac_num, treat_start):
    # EOT = End of Treatment
    frac_num = int(frac_num);
    weeks = frac_num//5;
    total_days = frac_num + 2*weeks;
    treat_days = np.tile([1,1,1,1,1,0,0], weeks+1)[:total_days].nonzero()[0]+treat_start;
    sim_resume_days = treat_days+10/(60*24); # ODE simulation resumes 10 minutes after RT
    return sim_resume_days[-1];

def get_first_repop(frac_num,treat_start):
    frac_num = int(frac_num);
    weeks = frac_num//5;
    total_days = frac_num + 2*weeks;
    treat_days = np.tile([1,1,1,1,1,0,0], weeks+1)[:total_days].nonzero()[0]+treat_start;
    sim_resume_days = treat_days+10/(60*24); # ODE simulation resumes 10 minutes after RT
    return sim_resume_days[0]

def get_schedule(frac_num,treat_start):
    frac_num = int(frac_num);
    weeks = frac_num//5;
    total_days = frac_num + 2*weeks;
    treat_days = np.tile([1,1,1,1,1,0,0], weeks+1)[:total_days].nonzero()[0]+treat_start;
    sim_resume_days = treat_days+10/(60*24); # ODE simulation resumes 10 minutes after RT
    return sim_resume_days,treat_days


# example calc:
#   frac_num = 5;25*(1+5/(8.5)) - (calc_EOT(frac_num,treat_start) - get_first_repop(frac_num,treat_start))/(0.17*3.9)
def get_EOT_index_pre(EOT,time_vec):
    return np.nonzero(np.equal(EOT, time_vec))[0][0];

def get_EOT_index(frac_num, treat_start, time_vec):
    return get_EOT_index_pre(calc_EOT(frac_num,treat_start), time_vec);

def parameter_setup(switch_vec, misc_pars):
    #weak_feedbackQ
    subSelectQ, use_muQ, compDosesQ, deathFdbkQ, c_dep_sQ, kimReprogQ, kimDeathValQ, kimICQ, model0Q, model1Q, model2Q = switch_vec;
    DT, post_therapy_end, pwr, h1, h2, l_vec, ss, a, b = misc_pars;
    
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
        v_suffix = ["2x30","2.67x15","3.4x10","5x5"];#["2x30","2.67x15","3.4x10","5x5"]; # 4 x 13
        v_suffix_save_A = "dose_comp";
        v_suffix_save = "dose_comp_"+str(post_therapy_end);#"dose_comp_0_reprog";
        # Doses = [1, 2, 40/15, 34/10, 5]; # dose fraction sizes in Gy
        # Frac = [60, 30, 15 ,10, 5]; 
        #Doses = np.arange(1,11).tolist();#, 34/10, 5]; # dose fraction sizes in Gy
        Doses = [2.0,6.0];#sorted(np.arange(1,21,dtype=float).tolist() + [40/15,34/10]);
        Frac = [30.0,7.0];#list(map(lambda d: float(np.floor(calc_BED_Frac(.17,.02,d))),Doses)); # 25*(1+(5)/8.5)
        # Doses = [2.0,2.0];#np.sort(np.append(np.arange(1,21),[40/15,34/10])).tolist(); # dose fraction sizes in Gy
        # Frac = [30.0,31.0];#list(map(lambda d:calc_BED_Frac(a,b,d),Doses));#,baseline_BED=25*(1+5/8.5)
        C = [0.0,5.196*10 ** (-pwr)];
        fg_subtitle = ['Dediff Rate:',str(C[0])];
    else:
        v_suffix = ["w/o reprog","w/ reprog"]; 
        v_suffix_save_A = "dediff_comp";
        v_suffix_save = "dediff_comp_"+str(post_therapy_end);
        Doses = [2,40/15,3.4,5]; # dose fraction sizes in Gy
        Frac = [30,15,10,5]; 
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
    cont_p_a = 0; 
    #chi = 0; 
    beta = 'tbd';
    xi1 = 'tbd';
    xi2 = 'tbd';
    mu_bar = 0; group_name = 'kim';
    if model0Q == True and model1Q == False and model2Q == False:
        group_name = 'm0';
    elif model0Q == False and model1Q == True and model2Q == False:
        cont_p_a = 10**4;
        group_name = 'm1';
    elif model0Q == False and model1Q == False and model2Q == True:
        cont_p_a = 10**4;
        chi = 10**4;
        mu_bar = 0.0143;
        group_name = 'm2';
    else:
        print("this ain't it chief")
    
    # will survivin be added?
    if use_muQ:
        mu_start = mu_bar;
    else:
        mu_start = 0;
    
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
        if type(ss) == int:
            rng = [ss];
        elif type(ss) == list or type(ss) == np.ndarray:
            rng = ss;
    else:
        rng = list(range(len(l_vec)));
    # no feedback, only div feedback, only prob feedback (weak, strong), 
    # both div and prob feedback (weak, strong)
    params = [total_start_frac, d, Doses, Frac, C, cont_p_a, beta, mu_bar, hd, cont_c, rng, mu_start, xi1, xi2, 'filler'];#[total_start_frac, d, Doses, Frac, C, cont_p_a, chi, mu_bar, hd, cont_c, rng, mu_start];
    stgs = [hd_str_mod, v_suffix, v_suffix_save_A, v_suffix_save, fg_subtitle, reprog_suffix, hd_str, hd_suffix, group_name, c_fdbk_mod];
    return [params, stgs]

def data_saver():
    return