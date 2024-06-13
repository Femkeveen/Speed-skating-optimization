#!/usr/bin/env python
# coding: utf-8

# In[78]:


def odetmt_symbolic():
    import numpy as np 
    import sympy as sm
    import sympy.physics.mechanics as me
    me.init_vprinting(use_latex='mathjax')
    t = me.dynamicsymbols._t

######################## Describing the SSM model symbolically ################  

    #absolute coordination 
    x_b, y_b, z_b, phi_b, x_s, y_s, z_s, phi_s = me.dynamicsymbols('x_b, y_b, z_b, phi_b, x_s, y_s, z_s, phi_s')
    x = sm.Matrix([x_b, y_b, z_b, x_s, y_s, phi_s])

    #generalized coordinates 
    theta_b, theta_s, u_s, v_s, w_s, u_b, v_b, w_b = me.dynamicsymbols('theta_b, theta_s, u_s, v_s, w_s, u_b, v_b, w_b')
    qd = sm.Matrix([u_b, v_b, w_s])     #generalized known coordinates
    qo = sm.Matrix([u_s, v_s, theta_s]) #generalized unknown coordinates 
    qi = sm.Matrix([qd, qo])            #generalized coordinates 

    #derivatives of generalized coordinates
    dqd = qd.diff(t) 
    dqo = qo.diff(t)
    dqi = sm.Matrix([dqd, dqo])

    #double derivatives 
    ddqo = dqo.diff(t)
    ddqd = dqd.diff(t)
    ddqi = dqi.diff(t)

    kk = sm.symbols('kk')          #sign deciding right (k = -1) or left (k = +1) skate
    m, m_b, m_s, I_s, M_s = sm.symbols('m, m_b, m_s, I_s, M_s') #mass, inertia and moment
    mu, g, k1, alpha = sm.symbols('mu, g, k1, alpha') #mu: ice friction constanct, m: mass of skater, g: gravitational acc

    constants = sm.Matrix([g, k1, kk, m, m_b, m_s, mu, alpha, I_s])

    #Reference frames 
    N = me.ReferenceFrame('N')     # Global reference frame 
    B = me.ReferenceFrame('B')     # Reference Frame of the body 
    S = me.ReferenceFrame('S')     # Reference Frame of the skate

    B.orient_axis(N, theta_b, N.z)    #orientation of the Body frame relative to the global, by heading of the body (theta_b)
    S.orient_axis(N, kk*theta_s, N.z) #orientation of the Skate frame relative to the global, ccw = + 

    ############################## Positions ########################################

    #creating points for the body and the skate 
    O = me.Point('O')  #origin 
    body = me.Point('body')
    skate = me.Point('skate')

    #set the position of the COM of the body relative to the origin in global coordinates 
    body.set_pos(O, x_b*N.x + y_b*N.y + z_b*N.z) 
    pos_body_x = body.pos_from(O)
    #Set the position of the skate relative to the origin in global coordinates 
    skate.set_pos(O, x_s*N.x+ y_s *N.y + z_s*N.z)
    pos_skate_x = skate.pos_from(O)
    #generalized coordinates of the skate, expressed in global coordinates. 
    us_expr = body.pos_from(skate).express(S).dot(S.x)
    vs_expr = body.pos_from(skate).express(S).dot(S.y)
    ws_exrp = body.pos_from(skate).express(S).dot(S.z)

    #position of body COM from origin in global reference frame in generalized coordinates
    body.set_pos(O, u_b*N.x + v_b*N.y + w_s*N.z)      #set position of the body in generalized coordinates
    pos_body_gen = body.pos_from(O)                   #position of the body from origin

    #positions of the skates from the body
    skate.set_pos(body,-u_s*S.y - v_s*kk*S.x-w_s*S.z) #define position of skate from body
    pos_skate_gen = skate.pos_from(body).express(N)   #position of skate from body expressed in global reference frame
    pos_skate_o_gen = skate.pos_from(O).express(N)    #position of skate form origin

    # # #Express x ([x_b, y_b, z_b, x_s, y_s, phi_s]) in generalized coordinates 
    x_expr = sm.Matrix([pos_body_gen.dot(N.x), pos_body_gen.dot(N.y), pos_body_gen.dot(N.z), (pos_body_gen+pos_skate_gen).dot(N.x), (pos_body_gen+pos_skate_gen).dot(N.y), kk*theta_s])
    x_dot = x.diff(t)
    x_expr_dot = x_expr.diff(t)

    ##Dict to translate from absolute to generalized coordinates (and derivatives)
    x_repl = dict(zip(x, x_expr))
    xd_repl = dict(zip(x_dot, x_expr_dot))
    x_repl, xd_repl


    ############################### Velocities ####################################
    # Velocities in absolute coordinates 
    O.set_vel(N,0) #set the velocity of the origin at 0
    vel_skate = skate.vel(N).express(S) #calculate the velocity of the skate and express it in reference frame of the skate

    lateral_vel_skate = vel_skate.dot(S.x) #calculate the velocity of the skate in the direction lateral to the heading of the skate
    lateral_vel_skate

    # # Translation to generalized coordinates
    later_vel_skate_gen = lateral_vel_skate.xreplace(xd_repl) #Express in generalized coordinates
    later_vel_skate_gen = sm.Matrix([later_vel_skate_gen])    #Translate to Matrix form
    later_vel_skate_gen #this is the non holonomic constraint 

    ############################# Transformation functions #########################

    T = x_expr_dot.jacobian(dqi) #Transformation matrix 

    Mij = sm.diag(m_b, m_b, m_b, m_s, m_s, I_s)  #Mass matrix 

    #Calculating the reduced mass matrix 
    Mredtot = sm.trigsimp(T.transpose()*Mij*T)  #TMT

    #Convective accelerations g 
    gcon  = x_expr_dot.jacobian(qi)*dqi

    # Jacobian of the non holonomic constraint (later_vel_skate) of the no-velocity condition in the lateral direction of the skate. 
    Cktot = later_vel_skate_gen.jacobian(dqi) #M_hn

    diff_constraint = later_vel_skate_gen.diff(t)  
    constraint_bias = diff_constraint.xreplace({qddr : 0 for qddr in ddqi}) #From https://moorepants.github.io/learn-multibody-dynamics/lagrange.html
    Cki = constraint_bias


#     ####################################### Forces ################################

    def friction(m, alpha, k1, mu, skate, x_b, y_b): 
        mb = (1-alpha)*m        #mass of the body 
        ms = alpha*m 

        vb = sm.sqrt(x_b.diff(t)**2 + y_b.diff(t)**2) #relative velocity of the speedskater to the air 
        F_b_friction = k1*vb**2

        F_N = m*g                 #approximation of normal force on ice 
        F_skate = F_N * mu        #ice friction on the skate 

        return F_b_friction, F_skate#, theta_b

    Fb, Fs = friction(m, alpha, k1, mu, skate, x_b, y_b)
    Fb = Fb.xreplace(x_repl)

    M_s = sm.symbols('M_s') # Moment around the skate

    #friction forces 
    F_fb = Fb*-B.y #friction force on body
    F_fs = Fs*-S.y #friction force on skate from ice 

    #gravitational forces on body
    Fz_b = m_b*g*-B.z 

    # Moment on the skate 
    M = kk * M_s

    #Force vector 
    fi = sm.Matrix([(F_fb+Fz_b).dot(N.x),   #total forces on body in global x dir 
                   (F_fb+Fz_b).dot(N.y),    #total forces on body in global y dir
                   (F_fb+Fz_b).dot(N.z),    #total forces on body in global z dir
                   -kk*F_fs.dot(N.x),       #total forces on skate in global x dir
                   F_fs.dot(N.y),           #total forces on skate in global y dir
                   M                        #Moment of skate around z axis. 
                  ])

    fi = fi.xreplace(xd_repl)               #express in generalized coordinates

    #Reduced force vector 
    fredtot = (T.transpose()*(fi-Mij*gcon))

    ############################## Solving equations ##############################
    #Calculating the solutions 
    Ltot = sm.Matrix([[Mredtot, Cktot.transpose()], [Cktot, sm.Matrix([0])]], )
    Ltot
    Rtot = -sm.Matrix([[fredtot], [-Cki]])
    Rtot

    #Reorganized in terms of known (q0) and unkown coodinates (qd) 
    locd = [0, 1, 2]
    loco = [3, 4, 5]
    Mdd = Mredtot[locd, :][:, locd]
    Mdo = Mredtot[locd, :][:, loco]
    Mod = Mredtot[loco, :][:, locd]
    Moo = Mredtot[loco, :][:, loco]
    Mdd, Mdo, Mod, Moo

    Ckd = Cktot[:, locd]
    Cko = Cktot[:, loco]
    Ckd, Cko

    Fd = fredtot[locd, :]
    Fo = fredtot[loco, :]

    M1 = sm.Matrix([[Mdd, Ckd.transpose()], [Ckd, sm.Matrix([0])]])
    M2 = sm.Matrix([[Fd - Mdo*ddqo], [-Cki - Cko*ddqo]])
    M1, M2

    sol1 = M1.LUsolve(M2) #solve [ddq, lambda]
    sol1

    ddxb = sol1[0]
    ddyb = sol1[1]
    ddzb = sol1[2]
    lamda = sol1[3]

    M3 = Mod.row_join(Moo).row_join(Cko.transpose()) # [Mod Moo Cko']
    M4 = sm.Matrix([[ddqd], [ddqo], [lamda]]) # [ddqd;ddqo;sol1(4)], expr 19
    sol2 = M3*M4 #[Mod Moo Cko']*[ddqd;ddqo;sol1(4)]=[Fo], Fo are forces on the skate
    sol2 #[Fd, Fo, -Ccon]


    ############################ Numerical evaluation of solutions ################
    ddus, ddvs, ddthetas = sm.symbols('ddus, ddvs, ddthetas')
    sol1_new = sol1.xreplace({ddqo[0]:ddus, ddqo[1]:ddvs, ddqo[2]: ddthetas})
#    sol1_new
# #     sol1_eval = sm.lambdify((g, k1, kk, m, m_b, m_s, mu, qi, dqi, theta_b, ddus, ddvs, ddthetas), sol1_new)
# #     sol1_num = sol1_eval(g_val, fric_coef, kk_val, m_skater, mb, ms, mu_val, qi_val, dqi_val, thetab_val, ddus_val, ddvs_val, ddthetas_val)
# #     print('sol1 = ', sol1_num)

    ddub, ddvb, ddws = sm.symbols('ddub, ddvb, ddws')
    ddxb_val, ddyb_val, ddzb_val = sol1[locd,:]
    
    ddqd_sol = sol1_new[locd,:]
#     ddqd_val[2] = ddws_val
#     ddqd_val, ddqo_val

    sol2_new = sol2.xreplace({ddqo[0]:ddus, ddqo[1]:ddvs, ddqo[2]: ddthetas, ddqd[0]:ddub, ddqd[1]:ddvb, ddqd[2]:ddws})
#     sol2_eval = sm.lambdify((I_s, g, k1, kk, m, m_b, m_s, mu, qi, dqi, theta_b, ddus, ddvs, ddthetas, ddub, ddvb, ddws), sol2_new)
#     sol2_num = sol2_eval(IS, g_val, fric_coef, kk_val, m_skater, mb, ms, mu_val, qi_val, dqi_val, thetab_val, ddus_val, ddvs_val, ddthetas_val, ddxb_val, ddyb_val, ddzb_val)
#     print('sol2 = ', sol2_num)

    ########################### Defining the outputs ##############################
    Fskate = sol2_new.reshape(3,1)
#     print('Fskate =', Fskate)
    lamda = sol1_new[3]
#     print('Lambda=', lamda_num)
#     def f_eval(function): 
#         f_eval = sm.lambdify((constants, theta_b, M_s, qi, dqi), function)
#         f_num = f_eval(constants_val, thetab_val, Ms, qi_val, dqi_val)
#         return f_num 
#     fredtot_num = f_eval(fredtot)
    Fi = fredtot
#     print('Fi =', Fi)
    Q = Fskate - Fi[loco,:]
#     print('Q=', Q)
    #dqd_val = np.array([dub_val, dvb_val, dws_val]).reshape(3,1)
    ydot = np.array([dqd, ddqd]).reshape(6,1)
#     ydot
#     print('ydot=', ydot)
#     dxs, dys = x_expr[3], x_expr[4]
  #  return ydot, 
    return sol1_new, sol2_new, ydot, lamda, Fskate, Fi, Q, Fs, Fb


# In[79]:


sol1, sol2, ydot_symb, lamda_symb, Fskate_symb, Fi_symb, Q_symb, Fs_symb, Fb_symb = odetmt_symbolic()

