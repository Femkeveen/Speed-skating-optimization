#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def rk4(b0,VAR,h,freqLPM,m_skater,alpha,fric_coef,mu,skate, VAR2): 
    from odetmt import odetmt
    import numpy as np
    #Calculate the forces with the TMT method 
    k1, Labdas, Forces, Q, Fi, F_wrijv, Fwrb, thetab = odetmt(b0, VAR[0], m_skater, alpha, fric_coef, mu, skate, freqLPM)
    #Update b0, 2nd order Runge Kutta? 
    b0new = b0 + (h / 2) * k1.T * (1 / freqLPM)

    #Calculate k2 
    k2, *_ = odetmt(b0new, VAR[1], m_skater, alpha, fric_coef, mu, skate, freqLPM)
    b0new = b0 + (h / 2) * k2.T * (1 / freqLPM)
    
    #Calculate k3 
    k3, *_ = odetmt(b0new, VAR[1], m_skater, alpha, fric_coef, mu, skate, freqLPM)
    b0new = b0 + h * k3.T * (1 / freqLPM)
    
    #Calculate k4
    k4, *_ = odetmt(b0new, VAR[2], m_skater, alpha, fric_coef, mu, skate, freqLPM)
    
    #update state variables bo = [xb, yb, zb,dxb, dyb, dzb]
    y = b0 + (h / 6) * (k1.T + 2 * k2.T + 2 * k3.T + k4.T) * (1 / freqLPM)
    y = y.reshape(6)
    #Calculate average slope K E
    K = (k1.T + 2 * k2.T + 2 * k3.T + k4.T) / 6 
    K = K.reshape(6)

    #Renaming variables 
    XB1, YB1, ZB1, dXB1, dYB1, dZB1 = y[0], y[1], y[2], y[3], y[4], y[5]   
    ddXB1, ddYB1, ddZB1 = K[3], K[4], K[5]
    
    mm = 2 #Taking the values at the end of the interval (3th row)
    # Determine variables based on skate side
    if skate == 'LS':
        US1, VS1, WS1, THETA_S1 = VAR[mm, 0], VAR[mm, 1], VAR[mm, 2], VAR[mm, 3]
        dUS1, dVS1, dWS1, dTHETA_S1 = VAR[mm, 8:12]
        ddUS1, ddVS1, ddWS1, ddTHETA_S1 = VAR[mm, 16:20]
    else:
        US1, VS1, WS1, THETA_S1 = VAR[mm, 4], VAR[mm, 5], VAR[mm, 6], VAR[mm, 7]
        dUS1, dVS1, dWS1, dTHETA_S1 = VAR[mm, 12:16]
        ddUS1, ddVS1, ddWS1, ddTHETA_S1 = VAR[mm, 20:24]
    
    #Calculate EPS and D based on constraints 
    #EPS is the non-holonomic constraint of no movement in the direction of the skate (see later_vel_skate_gen in odetmt). For kk = 1 and kk =-1
    # checked with LS = later_vel_skate_gen.xreplace({ kk : 1}), RS = later_vel_skate_gen.xreplace({kk:-1})
    if skate == 'LS':
        EPS = -np.sin(THETA_S1) * (-dUS1 * np.cos(THETA_S1) + dYB1 - dVS1 * np.sin(THETA_S1) - dTHETA_S1 * VS1 * np.cos(THETA_S1) + dTHETA_S1 * US1 * np.sin(THETA_S1)) - np.cos(THETA_S1) * (dXB1 - dVS1 * np.cos(THETA_S1) + dUS1 * np.sin(THETA_S1) + dTHETA_S1 * US1 * np.cos(THETA_S1) + dTHETA_S1 * VS1 * np.sin(THETA_S1))
        D = np.array([-np.cos(THETA_S1), -np.sin(THETA_S1), 0])
    else:
        EPS = -np.sin(THETA_S1) * (-dUS1 * np.cos(THETA_S1) + dYB1 - dVS1 * np.sin(THETA_S1) - dTHETA_S1 * VS1 * np.cos(THETA_S1) + dTHETA_S1 * US1 * np.sin(THETA_S1)) + np.cos(THETA_S1) * (-dUS1 * np.sin(THETA_S1) + dVS1 * np.cos(THETA_S1) + dXB1 - dTHETA_S1 * US1 * np.cos(THETA_S1) - dTHETA_S1 * VS1 * np.sin(THETA_S1))
        D = np.array([np.cos(THETA_S1), -np.sin(THETA_S1), 0])
    #D is the direction vector in which the non holonomic constraint force is applied 
    #(-S.x for right skate, S.x for left skate). And then express in global coordinates S.x.express(N)
    
    # Compute deltaX using Moore-Penrose pseudo-inverse, to ensure the new state adheres to the constraints
    deltaX = np.dot(D, np.linalg.pinv(np.outer(D, D))).dot(-EPS)

    # Update y with deltaX
    y += np.array([0, 0, 0, deltaX[0], deltaX[1], deltaX[2]])    

    #update state variables to new time step 
    ydotu, *_ = odetmt(y, VAR2, m_skater, alpha, fric_coef, mu, skate, freqLPM)

        
    return y,Labdas,Forces, Q, Fi, F_wrijv, Fwrb, ydotu, thetab 


