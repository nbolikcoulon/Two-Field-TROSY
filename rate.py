import numpy as np
from scipy import constants


import sys

#Constants for 3-fluorotyrosin



GAMMA_F = 2.51815e8
GAMMA_C = 6.72828e7
GAMMA_H = 2.67522e8

R_CF = 1.338e-10
R_FH = 2.61825e-10
R_CH = 2.14178e-10
R_HH = 4.27297e-10



MU4PI = constants.mu_0 / (4.0 * np.pi)

#PHI_CSA_DD_F = np.deg2rad(17.0)
#PHI_CSA_DD_C = np.deg2rad(16.0)

RATIO = GAMMA_C / GAMMA_F

#P2_CSA_DD_F = 0.5 * (3.0 * np.cos(PHI_CSA_DD_F) ** 2 - 1.0)
#P2_CSA_DD_C = 0.5 * (3.0 * np.cos(PHI_CSA_DD_C) ** 2 - 1.0)

AD_CF = -MU4PI * constants.hbar * GAMMA_F * GAMMA_C / (R_CF**3)
AD_CF_2 = AD_CF ** 2

AD_FH = -MU4PI * constants.hbar * GAMMA_F * GAMMA_H / (R_FH**3)
AD_FH_2 = AD_FH ** 2

AD_CH = -MU4PI * constants.hbar * GAMMA_C * GAMMA_H / (R_CH**3)
AD_CH_2 = AD_CH ** 2


def compute_spectral_density_Auto(w, s2, tauc, tauf):
    taui = tauc * tauf / (tauc + tauf)

    jw = 1.0/5.0 * (s2 * tauc / (1.0 + (w * tauc) ** 2) +
                      (1.0 - s2) * taui / (1.0 + (w * taui) ** 2))

    return jw





def compute_spectral_density_CSAcross(w, s2, tauc, tauf):
    taui = tauc * tauf / (tauc + tauf)

    jw = -1.0/2.0 * 1.0/5.0 * (s2 * tauc / (1.0 + (w * tauc) ** 2) +
                               (1.0 - s2) * taui / (1.0 + (w * taui) ** 2))

    return jw

def compute_spectral_density_CSAInterCross(w, s2, tauc, tauf, vecF, vecC):
    taui = tauc * tauf / (tauc + tauf)

    NormF = np.sqrt(vecF[0]**2 + vecF[1]**2 + vecF[2]**2)
    NormC = np.sqrt(vecC[0]**2 + vecC[1]**2 + vecC[2]**2)
    ScalarProd = vecF[0]*vecC[0] + vecF[1]*vecC[1] + vecF[2]*vecC[2]

    CosAngle = ScalarProd/(NormF*NormC)
    pS2 = (3.0 * CosAngle**2 - 1.0)/2.0

    jw = pS2 * 1.0/5.0 * (s2 * tauc / (1.0 + (w * tauc) ** 2) +
                               (1.0 - s2) * taui / (1.0 + (w * taui) ** 2))

    return jw


def compute_spectral_density_CrossCHF(w, s2, tauc, tauf):
    taui = tauc * tauf / (tauc + tauf)
    CosAngle = np.cos(0.546019)
    pS2 = (3.0 * CosAngle**2 - 1.0)/2.0

    jw = pS2 * 1.0/5.0 * (s2 * tauc / (1.0 + (w * tauc) ** 2) +
                               (1.0 - s2) * taui / (1.0 + (w * taui) ** 2))

    return jw



#def compute_rates(b0, s2, tauc, tauf, CSA_Val, vecFl, vecCl, Ftrosy, vecFp = [1, 0, 0], vecCp = [1, 0, 0], lambda_f_ext=0.0, rho_f_ext=0.0, rH = 0.0):
def compute_rates(b0, s2, tauc, tauf, CSA_Val, vecFl, vecCl, vecFp = [1, 0, 0], vecCp = [1, 0, 0], lambda_f_ext=0.0, rho_f_ext=0.0, rH = 0.0):
    rates = dict()
    
    
    CSA_F_long, CSA_F_perp, CSA_C_long, CSA_C_perp = CSA_Val[0], CSA_Val[1], CSA_Val[2], CSA_Val[3]
    
#    print(CSA_F_long, CSA_F_perp, CSA_C_long, CSA_C_perp)
#    print(vecFl)
#    print(vecCl)
#    sys.exit()
    
    

    w0_f = b0 * GAMMA_F
    w0_c = b0 * GAMMA_C
    w0_h = b0 * GAMMA_H

    ac_f_long = np.sqrt(2.0/3.0) * CSA_F_long * GAMMA_F * b0
    ac_c_long = np.sqrt(2.0/3.0) * CSA_C_long * GAMMA_C * b0
    ac_f_perp = np.sqrt(2.0/3.0) * CSA_F_perp * GAMMA_F * b0
    ac_c_perp = np.sqrt(2.0/3.0) * CSA_C_perp * GAMMA_C * b0
    

    omega_list = np.array([0.0, w0_c, w0_f - w0_c, w0_f, w0_f + w0_c])
    omega_list_FH = np.array([0.0, w0_h, w0_f - w0_h, w0_f, w0_f + w0_h])
    omega_list_CH = np.array([0.0, w0_h, w0_c - w0_h, w0_c, w0_c + w0_h])
    omega_list_CHF = np.array([0.0, w0_h])
    j0A, jwcA, jwmA, jwfA, jwpA = compute_spectral_density_Auto(omega_list, s2, tauc, tauf)
    j0AFH, jwHAFH, jwmAFH, jwfAFH, jwpAFH = compute_spectral_density_Auto(omega_list_FH, s2, tauc, tauf)
    j0ACH, jwHACH, jwmACH, jwcACH, jwpACH = compute_spectral_density_Auto(omega_list_CH, s2, tauc, tauf)
    j0CCH, jwHCCH = compute_spectral_density_CrossCHF(omega_list_CHF, s2, tauc, tauf)
    j0C, jwcC, jwmC, jwfC, jwpC = compute_spectral_density_CSAcross(omega_list, s2, tauc, tauf)
    j0FlCl, jwcFlCl, jwmFlCl, jwfFlCl, jwpFlCl = compute_spectral_density_CSAInterCross(omega_list, s2, tauc, tauf, vecFl, vecCl)
    j0FpCp, jwcFpCp, jwmFpCp, jwfFpCp, jwpFpCp = compute_spectral_density_CSAInterCross(omega_list, s2, tauc, tauf, vecFp, vecCp)
    j0FpCl, jwcFpCl, jwmFpCl, jwfFpCl, jwpFpCl = compute_spectral_density_CSAInterCross(omega_list, s2, tauc, tauf, vecFp, vecCl)
    j0FlCp, jwcFlCp, jwmFlCp, jwfFlCp, jwpFlCp = compute_spectral_density_CSAInterCross(omega_list, s2, tauc, tauf, vecFl, vecCp)
    j0FlCF, jwcFlCF, jwmFlCF, jwfFlCF, jwpFlCF = compute_spectral_density_CSAInterCross(omega_list, s2, tauc, tauf, vecFl, [0, 0, 1])
    j0ClCF, jwcClCF, jwmClCF, jwfClCF, jwpClCF = compute_spectral_density_CSAInterCross(omega_list, s2, tauc, tauf, vecCl, [0, 0, 1])
    
#    rates['AntiTrosyF'] = 1.0/12.0 * (
#            3.0 * AD_CF_2 * (4.0*j0A + 3.0*jwcA + jwmA + 3.0*jwfA + 6.0*jwpA) 
##            + AD_CH_2 * rH * (18.0*jwcACH + 6.0*jwmACH + 36.0*jwpACH)
##            + 6.0 * AD_FH_2 * (4.0*j0AFH + 3.0*jwfAFH + jwmAFH + 6.0*jwHAFH + 6.0*jwpAFH)
#            + (ac_c_long**2 + ac_c_perp**2) * 6.0 * jwcA + (ac_f_long**2 + ac_f_perp**2) * 2.0 * (4.0*j0A + 3.0*jwfA)
#            + ac_c_long*ac_c_perp * 12.0 * jwcC + ac_f_long*ac_f_perp * 4.0 * (4.0*j0C + 3.0*jwfC)
#            + Ftrosy * (np.sqrt(6.0) * AD_CF * 2.0 * ac_f_long * (4.0*j0FlCF + 3.0*jwfFlCF) + np.sqrt(6.0) * AD_CF * 2.0 * ac_f_perp * (4.0*j0C + 3.0*jwfC))
#            )
#
#    rates['TrosyF'] =1.0/12.0 * (
#            3.0 * AD_CF_2 * (4.0*j0A + 3.0*jwcA + jwmA + 3.0*jwfA + 6.0*jwpA) 
##            + AD_CH_2 * rH * (18.0*jwcACH + 6.0*jwmACH + 36.0*jwpACH)
##            + 6.0 * AD_FH_2 * (4.0*j0AFH + 3.0*jwfAFH + jwmAFH + 6.0*jwHAFH + 6.0*jwpAFH)
#            + (ac_c_long**2 + ac_c_perp**2) * 6.0 * jwcA + (ac_f_long**2 + ac_f_perp**2) * 2.0 * (4.0*j0A + 3.0*jwfA)
#            + ac_c_long*ac_c_perp * 12.0 * jwcC + (ac_f_long*ac_f_perp) * 4.0 * (4.0*j0C + 3.0*jwfC)
#            - Ftrosy * (np.sqrt(6.0) * AD_CF * 2.0 * ac_f_long * (4.0*j0FlCF + 3.0*jwfFlCF) + np.sqrt(6.0) * AD_CF * 2.0 * ac_f_perp * (4.0*j0C + 3.0*jwfC))
#            )
#
#    
#    rates['TrosyC'] = 1.0/12.0 * (
#            3.0*AD_CF_2 * (4.0*j0A + 3.0*jwcA + jwmA + 3.0*jwfA + 6.0*jwpA)
##            + AD_CH_2 * rH * 6.0 * (4.0*j0ACH + 3.0*jwcACH + jwmACH + 6.0*jwHACH + 6.0*jwpACH)
##            + AD_FH_2 * rH * (18.0*jwfAFH + 6.0*jwmAFH + 36.0*jwpAFH)
#            - np.sqrt(6.0) * AD_CF * 2.0 * ac_c_long * (    4.0*j0ClCF + 3.0*jwcClCF) - np.sqrt(6.0) * AD_CF * 2.0 * ac_c_perp * (4.0*j0C + 3.0*jwcC)
#            + 2.0 * (ac_c_long**2 + ac_c_perp**2) * (4.0*j0A + 3.0*jwcA) + 6.0 * (ac_f_long**2 + ac_f_perp**2) * jwfA
#            + ac_c_long*ac_c_perp * 4.0 * (4.0*j0C + 3.0*jwcC) + 12.0*ac_f_long*ac_f_perp * jwfC
#            )
    

    rates['R2f'] = 1.0/12.0 * ( 
        3.0 * AD_CF_2 * (4.0 * j0A + 6.0 * jwcA + jwmA + 3.0 * jwfA + 6.0 * jwpA)
#        + AD_FH_2 * rH * (24.0*j0AFH + 18.0*jwfAFH + 6.0*jwmAFH + 36.0*jwHAFH + 36.0*jwpAFH)
        + (ac_f_long**2 + ac_f_perp**2) * (8.0 * j0A + 6.0 * jwfA)
        + 2.0 * ac_f_long * ac_f_perp * (8.0 * j0C + 6.0 * jwfC)
    )
    

    rates['R2c'] = 1.0/12.0 * (
        3.0 * AD_CF_2 * (4.0 * j0A + 6.0 * jwfA + jwmA + 3.0 * jwcA + 6.0 * jwpA)
#        + AD_CH_2 * rH * (24.0*j0ACH + 18.0*jwcACH + 6.0*jwmACH + 36.0*jwHACH + 36.0*jwpACH)
        + (ac_c_long**2 + ac_c_perp**2) * (8.0 * j0A + 6.0 * jwcA)
        + 2.0 * ac_c_long * ac_c_perp * (8.0 * j0C + 6.0 * jwcC)
    )
    

    rates['rho_f'] = ( 
        AD_CF_2/2.0 * (jwmA + 3.0 * jwfA + 6.0 * jwpA)
        + AD_FH_2 * rH * (jwmAFH + 3.0*jwfAFH + 6.0*jwpAFH)
        + (ac_f_long**2 + ac_f_perp**2) * jwfA +
        2 * ac_f_long * ac_f_perp * jwfC
    )
    

    rates['rho_c'] = (
        AD_CF_2/2.0 * (jwmA + 3.0 * jwcA + 6.0 * jwpA)
        + AD_CH_2 * rH * (jwmACH + 3.0*jwcACH + 6.0*jwpACH)
        + (ac_c_long**2 + ac_c_perp**2) * jwcA +
        2 * ac_c_long * ac_c_perp * jwcC
    )


    rates['rhoa_f'] = 1.0/12.0 * (
        3.0 * AD_CF_2 * (4.0 * j0A + jwmA + 3.0 * jwfA + 6.0 * jwpA) +
        12.0* AD_FH_2 * rH * (2.0*j0AFH + 3.0/2.0*jwfAFH + 1.0/2.0*jwmAFH + 18.0*jwHAFH + 18.0*jwpAFH) +
        12.0* AD_CH_2 * rH * (3.0*jwcACH + jwmACH + 6.0*jwpACH) + 
        2.0 * (6.0 * (ac_c_long**2 + ac_c_perp**2) * jwcA + (ac_f_long**2 + ac_f_perp**2) * (4.0 * j0A + 3.0 * jwfA) +
        12.0 * ac_c_long * ac_c_perp * jwcC + ac_f_long * ac_f_perp * (8.0 * j0C + 6.0 * jwfC))
    )
    

    rates['rhoa_c'] = 1.0/12.0 * (
        3.0 * AD_CF_2 * (4.0 * j0A + jwmA + 3.0 * jwcA + 6.0 * jwpA) +
        12.0* AD_CH_2 * rH * (3.0/2.0*jwcACH + 2.0*j0ACH + 1.0/2.0*jwmACH + 3.0*jwHACH + 3.0*jwpACH) +
        12.0* AD_FH_2 * rH * (jwmAFH + 6.0*jwpAFH + 3.0*jwfAFH) + 
        2.0 * (6.0 * (ac_f_long**2 + ac_f_perp**2) * jwfA + (ac_c_long**2 + ac_c_perp**2) * (4.0 * j0A + 3.0 * jwcA) +
        12.0 * ac_f_long * ac_f_perp * jwfC + ac_c_long * ac_c_perp * (8.0 * j0C + 6.0 * jwcC))
    )
    
    

    rates['lambda_mq'] = 1.0/12.0 * (
        3.0 * AD_CF_2 * (3.0 * jwcA + jwmA + 3.0 * jwfA + 6.0 * jwpA) +
        AD_FH_2 * rH * (24.0*j0AFH + 18.0*jwfAFH + 6.0*jwmAFH + 36.0*jwHAFH +36.0*jwpAFH) + 
        6.0 * AD_CH_2 * rH * (4.0*j0ACH + 3.0*jwcACH + jwmACH + 6.0*jwHACH + 6.0*jwpACH) +
        (ac_c_long**2 + ac_c_perp**2) * (8.0 * j0A + 6.0 * jwcA) + (ac_f_long**2 + ac_f_perp**2) * (8.0 * j0A + 6.0 * jwfA) +
        ac_c_long*ac_c_perp * (16.0 * j0C + 12.0 * jwcC) + ac_f_long*ac_f_perp * (16.0 * j0C + 12.0 * jwfC)
    )
    

    rates['rho_2sp_cf'] = (
        3.0/2.0 * AD_CF_2 * (jwfA + jwcA) +
        AD_CH_2 * rH * (jwmACH + 6.0*jwpACH) + AD_FH_2 * rH * (jwmAFH + 6.0*jwpAFH) + 
        (ac_f_long**2 + ac_f_perp**2) *jwfA + (ac_c_long**2 + ac_c_perp**2) *jwcA + 
        2.0 * ac_f_long*ac_f_perp * jwfC + 2.0 * ac_c_long*ac_c_perp * jwcC
    )
    

    rates['sigma'] = 1.0/2.0 * AD_CF_2 * (-jwmA + 6.0 * jwpA)
    

    rates['mu_mq'] = -1.0/12.0 * (
#            6.0 * AD_CH*AD_FH * rH * (9.0*jwHCCH + 6.0*j0CCH) + 
            3.0 * AD_CF_2 * (jwmA - 6.0 * jwpA) - 
            16.0 * (ac_f_long*ac_c_long*j0FlCl + ac_f_long*ac_c_perp*j0FlCp + ac_f_perp*ac_c_long*j0FpCl + ac_f_perp*ac_c_perp*j0FpCp)
    )

    rates['delta_f'] = 2.0 * np.sqrt(3.0/2.0) * AD_CF * (ac_f_perp * jwfC + ac_f_long * jwfFlCF)
            
    rates['delta_c'] = 2.0 * np.sqrt(3.0/2.0) * AD_CF * (ac_c_perp * jwcC + ac_c_long * jwcClCF)

    rates['eta_f'] = 1.0/np.sqrt(6) * AD_CF * (
        ac_f_long * (4.0 * j0FlCF + 3.0 * jwfFlCF) + ac_f_perp * (4.0 * j0C + 3.0 * jwfC)
    )
    
    rates['eta_c'] = 1.0/np.sqrt(6) * AD_CF * (
        ac_c_long * (4.0 * j0ClCF + 3.0 * jwcClCF) + ac_c_perp * (4.0 * j0C + 3.0 * jwcC)
    )

    rates['theta_c'] = 1e-10 * (b0 * GAMMA_C * rates['rho_c'] + b0 * GAMMA_F * rates['sigma'])
    rates['theta_f'] = 1e-10 * (b0 * GAMMA_F * rates['rho_f'] + rates['sigma'] * b0 * GAMMA_C)
    rates['theta_cf'] = 1e-10 * (rates['delta_c'] * b0 * GAMMA_C + b0 * GAMMA_F * rates['delta_f'])
    
    
    

    return rates


#