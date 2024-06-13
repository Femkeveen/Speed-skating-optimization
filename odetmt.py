#!/usr/bin/env python
# coding: utf-8

# In[170]:


def odetmt(b0, q, m_skater, alpha_val, fric_coef, mu_val, skate_val, freqMS):
    import numpy as np
    b0 = b0.reshape(6)
    m_skater = m_skater
    g_val = 9.81
    alpha_val = 0.1


    mb = (1-(alpha_val))*m_skater  #mass of body COM
    ms = (alpha_val)*m_skater      #mass of skate
    Ms = 0;
    Ib = 0;
    IS = 0

    xb_val, yb_val, zb_val = b0[0], b0[1], b0[2]
    ub_val, vb_val = b0[0], b0[1]  

    dxb_val, dyb_val, dzb_val = b0[3], b0[4], b0[5]
    dub_val, dvb_val  = b0[3], b0[4]

    thetab_val = -np.arcsin(dxb_val / np.sqrt(dxb_val ** 2 + dyb_val ** 2))

    if skate_val == 'LS': #Left skate 
        us_val, vs_val, ws_val, thetas_val = q[0], q[1], q[2], q[3]
        dus_val, dvs_val, dws_val,  dthetas_val = q[8], q[9], q[10], q[11]
        ddus_val, ddvs_val, ddws_val, ddthetas_val = q[16], q[17], q[18], q[19]
        kk = 1 

    elif skate_val == 'RS': #Right skate 
        us_val, vs_val, ws_val, thetas_val = q[4], q[5], q[6], q[7]
        dus_val, dvs_val, dws_val,  dthetas_val = q[12], q[13], q[14], q[15]
        ddus_val, ddvs_val, ddws_val, ddthetas_val  = q[20], q[21], q[22], q[23]
        kk = -1

    #Organizing the input values   
    qi_val = np.array([ub_val, vb_val, ws_val, us_val, vs_val, thetas_val])         #generalized coordinates
    dqi_val = np.array([dub_val, dvb_val, dws_val, dus_val, dvs_val, dthetas_val])  #derivatives of gen coord
    qd_val = np.array([ub_val, vb_val, ws_val]) #generalized known coordinates
    dqd_val = np.array([dub_val, dvb_val, dws_val])

    qo_val = np.array([us_val, vs_val, thetas_val]) 
    dqo_val = np.array([dus_val, dvs_val, dthetas_val])
    ddqo_val = np.array([ddus_val, ddvs_val, ddthetas_val])
    ddqo_val = np.array([ddus_val, ddvs_val, ddthetas_val])

    #dxs = -kk * vs_val * np.cos(kk * thetas_val) + ub_val + us_val * np.sin(kk * thetas_val)

    Mij = np.diag([mb, mb, mb, ms, ms, IS])
    thetas, dthetas = thetas_val, dthetas_val
    us, vs, dus, dvs  = us_val, vs_val, dus_val, dvs_val
    vb, dvb, ub, dub = vb_val, dvb_val, ub_val, dub_val
    
    #Transformation matrix 
    T =  np.array([
    [1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [1, 0, 0, np.sin(kk*thetas), -kk * np.cos(kk*thetas), kk * us * np.cos(kk*thetas) + kk**2 * vs * np.sin(kk*thetas)],
    [0, 1, 0, -np.cos(kk*thetas), -kk*np.sin(kk*thetas), kk*us*np.sin(kk*thetas) - kk**2*vs*np.cos(kk*thetas)],
    [0, 0, 0, 0, 0, kk]
    ])
    T

    #Convective acceleration
    gcon = np.array([
        0,
        0,
        0,
        dthetas * (dthetas * (kk * vs * np.cos(thetas) - kk * us * np.sin(thetas)) + dus * kk * np.cos(thetas) + dvs * kk * np.sin(thetas)) + 
        dthetas * dus * kk * np.cos(thetas) + dthetas * dvs * kk * np.sin(thetas),
        dthetas * (dthetas * (us * np.cos(thetas) + vs * np.sin(thetas)) - dvs * np.cos(thetas) + dus * np.sin(thetas)) - 
        dthetas * dvs * np.cos(thetas) + dthetas * dus * np.sin(thetas),
        0
    ])
    
    Cki = dthetas * (dvs * (-2 * np.cos(thetas) * np.sin(thetas) * kk**2 + 2 * np.cos(thetas) * np.sin(thetas)) -
                 dvb * np.cos(thetas) + 
                 dthetas * (np.cos(thetas) * (vs * np.cos(thetas) - us * np.sin(thetas)) -
                             np.sin(thetas) * (us * np.cos(thetas) + vs * np.sin(thetas)) -
                             kk * np.cos(thetas) * (kk * vs * np.cos(thetas) - kk * us * np.sin(thetas)) +
                             kk * np.sin(thetas) * (kk * us * np.cos(thetas) + kk * vs * np.sin(thetas))) + \
     dus * (np.cos(thetas)**2 - np.sin(thetas)**2 - kk**2 * np.cos(thetas)**2 + kk**2 * np.sin(thetas)**2) + \
     dub * kk * np.sin(thetas)) - \
     dthetas * dus * (kk**2 * np.cos(thetas)**2 + np.sin(thetas)**2) + \
     dthetas * dvs * (-np.cos(thetas) * np.sin(thetas) * kk**2 + np.cos(thetas) * np.sin(thetas))



    Cktot = np.array([
        np.cos(kk*thetas), 
        np.sin(kk*thetas), 
        0, 
        0, 
        -kk*np.sin(kk*thetas)**2 - kk*np.cos(kk*thetas)**2,
        (kk**2 * vs * np.sin(kk * thetas) + kk*us*np.cos(kk*thetas))*np.cos(kk * thetas) + (-kk**2 * vs * np.cos(kk * thetas) + kk * us * np.sin(kk * thetas)) * np.sin(kk * thetas)

    ])

    #Friction forces 
    Fb = fric_coef*np.sqrt(dub**2+dvb**2)**2
    Fs = g_val*mu_val*m_skater

    #Force vector 
    fi = np.array([
        Fb*np.sin(thetab_val), 
        -Fb*np.cos(thetab_val), 
        -g_val*mb,
        -kk*Fs*np.sin(kk*thetas), 
        -Fs*np.cos(kk*thetas),
        Ms*kk
    ])

    #Reduced mass matrix and force vector 
    #Mredtot = np.dot(np.dot(T.transpose(), Mij), T)
    Mredtot = np.matmul(np.matmul(T.T, Mij), T)
    fredtot = np.dot(T.T, fi - np.dot(Mij, gcon))
    Fi = np.dot(T.T, fi - np.dot(Mij, gcon))

    #Reorganization in terms of known and unkown coordinates 
    locd = [0, 1, 2]
    loco = [3, 4, 5]
    Mdd = Mredtot[locd, :][:, locd]
    Mdo = Mredtot[locd, :][:, loco]
    Mod = Mredtot[loco, :][:, locd]
    Moo = Mredtot[loco, :][:, loco]
    Mdd, Mdo, Mod, Moo

    Ckd = Cktot[locd]
    Cko = Cktot[loco]
    Ckd, Cko

    Fd = fredtot[locd]
    Fo = fredtot[loco]
    Fd, Fo


    Ckd_row = np.array([[Ckd[0], Ckd[1], Ckd[2]]])
    Cko_row = np.array([[Cko[0], Cko[1], Cko[2]]])

    ########################## Solving the equations ##################################
    # [Mdd Mdo Cd.T; Ckd Cko.T 0] * [ddqd ddqo lambda]' = Fd
    
    #Left hand side matrix 
    M1 = np.block([[Mdd, Ckd_row.T],
                   [Ckd_row, np.zeros((1, 1))]])

    # # Construct the right-hand side vector
    #M2 = np.vstack((Fd - np.dot(Mdo, ddqo_val), -Cki - np.dot(Cko, ddqo_val)))
    e = Fd - np.dot(Mdo, ddqo_val) 
    M2 = np.array([
        [e[0], e[1], e[2], -Cki - np.dot(Cko, ddqo_val)]
    ])
    # # Solve the linear system
    sol1 = np.linalg.solve(M1, M2.transpose())

    # print("sol1:", sol1)``b
    sol1

    ddqd = sol1[locd]
    ddqd[2] = ddws_val
    
    lamda = sol1[3]
    
    # [Mod Moo Cko] * [ddqd ddqo lamda] = [Fo]

    M3 = np.hstack((Mod, Moo, Cko_row.T)) 
    M4 = np.hstack((ddqd.reshape(3), ddqo_val, sol1[3]))
    sol2 = np.dot(M3, M4)
    
    Fskate = sol2.reshape(3)
    Fi = Fi.reshape(6)
    
    Q = Fskate - Fi[loco]
    Q = Q.reshape(3)
    ydot = np.hstack((dqd_val, ddqd.reshape(3)))
    ydot = ydot.reshape(6)

    
    return ydot, lamda, Fskate, Q, Fi, Fs, Fb, thetab_val





