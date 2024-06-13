#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def gen_coord_meas(POS, VEL, T):
    import numpy as np

    THETA_LS, THETA_RS, VLS, ULS, WLS, VRS, URS, WRS, THETA_B = np.zeros(T), np.zeros(T), np.zeros(T), np.zeros(T), np.zeros(T), np.zeros(T), np.zeros(T), np.zeros(T), np.zeros(T) 

    for k in range(T):
            xb = POS[0, k]
            yb = POS[1, k]
            zb = POS[2, k]
            xls = POS[3, k]
            yls = POS[4, k]
            zls = POS[5, k]
            xrs = POS[6, k]
            yrs = POS[7, k]
            zrs = POS[8, k]

            dxb = VEL[0, k]
            dyb = VEL[1, k]
            dzb = VEL[2, k]
            dxls = VEL[3, k]
            dyls = VEL[4, k]
            dzls = VEL[5, k]
            dxrs = VEL[6, k]
            dyrs = VEL[7, k]
            dzrs = VEL[8, k]

            # Calculate THETA_B
            THETA_B[k] = -np.arcsin(dxb / np.sqrt(dyb**2 + dxb**2))

            # Calculate THETA_LS
            THETA_LS[k] = -np.arcsin(dxls / np.sqrt(dyls**2 + dxls**2))
            THETA_RS[k] = np.arcsin(dxrs / np.sqrt(dyrs**2 + dxrs**2))

            # Calculate VLS, ULS, WLS
            VLS[k] = np.sin(THETA_LS[k]) * (yb - yls) + np.cos(THETA_LS[k]) * (xb - xls)
            ULS[k] = np.cos(THETA_LS[k]) * (yb - yls) - np.sin(THETA_LS[k]) * (xb - xls)
            WLS[k] = zb - zls

            # Calculate VRS, URS, WRS
            VRS[k] = np.sin(THETA_RS[k]) * (yb - yrs) - np.cos(THETA_RS[k]) * (xb - xrs)
            URS[k] = np.cos(THETA_RS[k]) * (yb - yrs) + np.sin(THETA_RS[k]) * (xb - xrs)
            WRS[k] = zb - zrs

    return THETA_B, THETA_LS, THETA_RS, VLS, ULS, WLS, VRS, URS, WRS

