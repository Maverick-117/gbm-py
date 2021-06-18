# -*- coding: utf-8 -*-
"""
Created on Wed May 19 14:23:47 2021

@author: jhvo9
"""
from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np

def dX_dt(X, t, P, l, h, m_u, m_v, a_v, n):
    return np.array([(2*P/(1+l*X[1]**n)-1)*m_u/(1+h*X[1]**n)*X[0], 
                    2*(1-P/(1+l*X[1]**n))*m_u/(1+h*X[1]**n)*X[0] + X[1]*(m_v/(1+h*X[1]**n) - a_v)])


def dynamo(sim_pars, dX_pars):
    total_start_frac, F, resume_rad, a1, a2, b1, b2, c, D, treat_days = sim_pars;
    P, l, h, m_u, m_v, a_v, n = dX_pars;
    t = np.linspace(0, 100, 2000, endpoint=False)
    X0 = np.array([total_start_frac*F, total_start_frac*(1-F)])
    sol = integrate.odeint(dX_dt, X0, t, args=(P, l, h, m_u, m_v, a_v, n)).T
    u, v = sol;
    
    # start from last term of u,v
    for i in range(len(resume_rad)):
        u, v = u[-1], v[-1]
        # switch terms around
        u, v = u*np.exp(-a1*D-b1*D**2) + c*D*v, v*np.exp(-a2*D-b2*D**2) - c*D*v
        
        # define new time interval, 1000 dense per day
        t_var = np.linspace(resume_rad[i], treat_days[i+1], int(np.round(1000*(-resume_rad[i] + treat_days[i+1]))))
          
        # solve ODE
        current_sol = integrate.odeint(dX_dt, [u, v], t_var, args=(P, l, h, m_u, m_v, a_v, n)).T
        
        # update sol, t
        u, v = current_sol
        sol = np.hstack((sol, current_sol))
        t = np.concatenate((t, t_var))
    
    return sol, t
      
  