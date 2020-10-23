
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 15:09:01 2019

@author: nbc
"""
import numpy as np

import rate as R
import liouvillian as L

import sys

from matplotlib import pyplot as plt

#np.set_printoptions(precision=2)
#np.set_printoptions(suppress=True)


GAMMA_F = 2.51815e8
GAMMA_C = 6.72828e7

#Pulses
#HF
tau_90_F = 10e-6
tau_90_C = 12e-6
b1_F = (np.pi/2.0) * 1.0/(GAMMA_F * tau_90_F)
b1_C = (np.pi/2.0) * 1.0/(GAMMA_C * tau_90_C)

tau_180_F = 2.0*tau_90_F
tau_180_C = 2.0*tau_90_C

tau_240_F = 240.0/90.0*tau_90_F




#LF
tau_90_F_LF = 6e-6
tau_90_C_LF = 10e-6
b1_F_LF = (np.pi/2.0) * 1.0/(GAMMA_F * tau_90_F_LF)
b1_C_LF = (np.pi/2.0) * 1.0/(GAMMA_C * tau_90_C_LF)

tau_180_F_LF = 2.0*tau_90_F_LF
tau_180_C_LF = 2.0*tau_90_C_LF






JCF = -240.0 #in Hz



def SimulateDummyScans(ThermalHF, ThermalLF, s2, tauc, tauf, ChemShift_Val, CSA_val, vecFl, vecCl, vecFp, vecCp, tmax, incr, b0, b0LF, delays, FieldList, nH, d1):

    bUp, aUp, bDown, aDown, incrShuttle = delays
    
    #INEPT evolution time
    T = 1.0/(2.0*np.abs(JCF))
    
    t1max, t2max = tmax
    incr1, incr2 = incr
    t1Delays = np.arange(0.0, t1max+incr1, incr1)
    
    
    
    
######### HIGH FIELD PULSES
#Carbon 90-pulses propagators    
    P_90c_px = L.make_pulse_c(ThermalHF, w1_c=b1_C, phase_c=0.0, pw=tau_90_C)
    P_90c_py = L.make_pulse_c(ThermalHF, w1_c=b1_C, phase_c=1.0, pw=tau_90_C)
    
    P_90c_mxpy = L.make_pulse_c(ThermalHF, w1_c=b1_C, phase_c=1.0/2.0, pw=tau_90_C)
    
    
    
#Simulatious 180-pulses propagators
    P_180fc_px = L.make_simpulse_cf(ThermalHF, w1_f=b1_F, phase_f=1.0, pw_f=tau_180_F, w1_c=b1_C, phase_c=0.0, pw_c=tau_180_C)
    P_180f_px = L.make_pulse_f(ThermalHF, w1_f=b1_F, phase_f=1.0, pw=tau_180_F)
    
    
    
######### LOW FIELD PULSES
#Fluor pulses propagators
    P_90f_px_LF = L.make_pulse_f(ThermalLF, w1_f=b1_F_LF, phase_f=0.0, pw=tau_90_F_LF)
    P_90f_py_LF = L.make_pulse_f(ThermalLF, w1_f=b1_F_LF, phase_f=1.0, pw=tau_90_F_LF)
    
    P_180f_px_LF = L.make_pulse_f(ThermalLF, w1_f=b1_F_LF, phase_f=0.0, pw=tau_180_F_LF)

#Carbon 180 pulse
    P_180c_px_LF = L.make_pulse_c(ThermalLF, w1_c=b1_C_LF, phase_c=0.0, pw=tau_180_C_LF)
    
    
    
######### SPIN FREE EVOLUTION
#INEPT spin evolution
    P_S3E_HF = L.make_delay(ThermalHF, T/4.0 - max(tau_240_F, tau_180_C)/np.pi - tau_90_F*2.0/np.pi)
    
#Equilibration time
    P_Eq = L.make_delay(ThermalHF, 1e9)
    
    
#d1 evolve
    P_d1 = L.make_delay(ThermalHF, t2max+d1)
    
    
    
#Waiting time
    BeforeUp = L.make_delay(ThermalHF, bUp)
    AfterUp = L.make_delay(ThermalLF, aUp)
    
    BeforeDown = L.make_delay(ThermalLF, bDown)
    AfterDown = L.make_delay(ThermalHF, aDown)
    
#Shuttling
    #Selector
    SelectionShuttlingUp = np.zeros((16, 16), dtype=np.complex128)
    SelectionShuttlingUp[0, 0] = 1.0
    SelectionShuttlingUp[3, 3] = 1.0
    
    SelectionShuttlingDown = np.zeros((16, 16), dtype=np.complex128)
    SelectionShuttlingDown[0, 0] = 1.0
    SelectionShuttlingDown[3, 3] = 1.0
    SelectionShuttlingDown[6, 6] = 1.0
    SelectionShuttlingDown[15, 15] = 1.0
    
    
    
    Prop_Up = [[] for F in FieldList]
    Prop_Down = [[] for F in FieldList]
    for F in range(len(FieldList)):
        NewRates = R.compute_rates(FieldList[F], s2, tauc, tauf, CSA_val, vecFl, vecCl, vecFp = [1, 0, 0], vecCp = [1, 0, 0], lambda_f_ext=0.0, rho_f_ext=0.0, rH = nH/2.0)
        
        NewL = L.compute_free_precess_liouv(ChemShift_Val[0], ChemShift_Val[1], FieldList[F], NewRates, JCF)
        Thermal_down = L.ThermalCorr(NewL, NewRates)
        
        NewL_up =  L.ThermalCorr(NewL, NewRates)
        
        Prop_Up[F] = L.make_delay(NewL_up, incrShuttle)
        Prop_Down[F] = L.make_delay(Thermal_down, incrShuttle)
        
    P_up = Prop_Up[0]
    P_down = Prop_Down[0]
    for F in range(len(FieldList)):
        P_up = Prop_Up[F] @ P_up
        P_down = P_down @ Prop_Down[F]
        
        
        
    
#t1 LF evolution
    P_t1 = []

    ta = T/2.0 - tau_180_C_LF/np.pi
    tb = 0.0 - tau_180_F_LF/np.pi - tau_180_C_LF/np.pi
    tc = T/2.0 - tau_180_F_LF/np.pi
    
    
    incr_a = incr1/2.
    incr_b = incr1/2. - T/(2.*(len(t1Delays)-1.))
    incr_c = T/(2.*(len(t1Delays)-1.))
    
    
    for t in range(len(t1Delays)):
        
        La = L.make_delay(ThermalLF, ta)
        Lb = L.make_delay(ThermalLF, tb)
        Lc = L.make_delay(ThermalLF, tc)
        
        P_t1.append(Lc @ P_180f_px_LF @ Lb @ P_180c_px_LF @ La)
        
        ta += incr_a
        tb += incr_b
        tc -= incr_c

        
######### PHASE CYCLING
    Phi1_x = P_90f_px_LF
    
    
    Phi2 = P_90f_py_LF
    Phi3 = P_90c_mxpy
    
    
    
######### SIMULATE SEQUENCE
    M_LF = SelectionShuttlingUp @ AfterUp @ P_up @ BeforeUp
    
#Low field evolution
    M_PostT1_x = Phi2 @ P_t1[0] @ Phi1_x @ M_LF

    PDOWN = P_180f_px @ AfterDown @ P_down @ SelectionShuttlingDown @ BeforeDown

    
#High field evolution
    Mdetect = P_90c_py @ Phi3 @ P_S3E_HF @ P_180fc_px @ P_S3E_HF @ P_90c_px @ PDOWN
    
    Mtot = P_d1 @ Mdetect @ M_PostT1_x
    
    
    Receiver_Longitunial = [np.zeros((1, 16), dtype=np.complex128) for i in range(3)]
    Receiver_Longitunial[0][0, 3] = 1.0     #Fz
    

    Prop = P_Eq @ L.set_fzcz_eq()
    for n in range(10):
        M = Mtot @ Prop
        Prop = M
        
        
    return Prop






def CdetecTROSY_S3E_ST2PT(s2, tauc, tauf, ChemShift_Val, CSA_val, vecFl, vecCl, vecFp, vecCp, tmax, incr, b0, b0LF, delays, FieldList, nH, d1):

    bUp, aUp, bDown, aDown, incrShuttle = delays
    
    #INEPT evolution time
    T = 1.0/(2.0*np.abs(JCF))
    
    t1max, t2max = tmax
    incr1, incr2 = incr
    t1Delays = np.arange(0.0, t1max+incr1, incr1)
    t2Delays = np.arange(0.0, t2max+incr2, incr2)
    
    
    
    
    ratesHF = R.compute_rates(b0, s2, tauc, tauf, CSA_val, vecFl, vecCl, vecFp = [1, 0, 0], vecCp = [1, 0, 0], lambda_f_ext=0.0, rho_f_ext=0.0, rH = nH/2.0)
    ratesLF = R.compute_rates(b0LF, s2, tauc, tauf, CSA_val, vecFl, vecCl, vecFp = [1, 0, 0], vecCp = [1, 0, 0], lambda_f_ext=0.0, rho_f_ext=0.0, rH = nH/2.0)
    LiouvillianHF = L.compute_free_precess_liouv(ChemShift_Val[0], ChemShift_Val[1], b0, ratesHF, JCF)
    LiouvillianLF = L.compute_free_precess_liouv(ChemShift_Val[0], ChemShift_Val[1], b0LF, ratesLF, JCF)
    
    

    ThermalHF = L.ThermalCorr(LiouvillianHF, ratesHF)
    ThermalLF = L.ThermalCorr(LiouvillianLF, ratesLF)
    
######### HIGH FIELD PULSES
#Carbon 90-pulses propagators    
    P_90c_px = L.make_pulse_c(ThermalHF, w1_c=b1_C, phase_c=0.0, pw=tau_90_C)
    P_90c_py = L.make_pulse_c(ThermalHF, w1_c=b1_C, phase_c=1.0, pw=tau_90_C)
    
    P_90c_mxpy = L.make_pulse_c(ThermalHF, w1_c=b1_C, phase_c=1.0/2.0, pw=tau_90_C)
    P_90c_mxmy = L.make_pulse_c(ThermalHF, w1_c=b1_C, phase_c=5.0/2.0, pw=tau_90_C)
    
    
#Simulatious 180-pulses propagators
    P_180fc_px = L.make_simpulse_cf(ThermalHF, w1_f=b1_F, phase_f=1.0, pw_f=tau_180_F, w1_c=b1_C, phase_c=0.0, pw_c=tau_180_C)
    P_180f_px = L.make_pulse_f(ThermalHF, w1_f=b1_F, phase_f=1.0, pw=tau_180_F)
    
    
    
######### LOW FIELD PULSES
#Fluor pulses propagators
    P_90f_px_LF = L.make_pulse_f(ThermalLF, w1_f=b1_F_LF, phase_f=0.0, pw=tau_90_F_LF)
    P_90f_mx_LF = L.make_pulse_f(ThermalLF, w1_f=b1_F_LF, phase_f=2.0, pw=tau_90_F_LF)
    P_90f_py_LF = L.make_pulse_f(ThermalLF, w1_f=b1_F_LF, phase_f=1.0, pw=tau_90_F_LF)
    P_90f_my_LF = L.make_pulse_f(ThermalLF, w1_f=b1_F_LF, phase_f=3.0, pw=tau_90_F_LF)
    
    P_180f_px_LF = L.make_pulse_f(ThermalLF, w1_f=b1_F_LF, phase_f=0.0, pw=tau_180_F_LF)

#Carbon 180 pulse
    P_180c_px_LF = L.make_pulse_c(ThermalLF, w1_c=b1_C_LF, phase_c=0.0, pw=tau_180_C_LF)
    
    
    
######### SPIN FREE EVOLUTION
#INEPT spin evolution
    P_S3E_HF = L.make_delay(ThermalHF, T/4.0 - max(tau_240_F, tau_180_C)/np.pi - tau_90_F*2.0/np.pi)
    
#Equilibration time
#    P_Eq = L.make_delay(ThermalHF, 1e3)
    P_steady = SimulateDummyScans(ThermalHF, ThermalLF, s2, tauc, tauf, ChemShift_Val, CSA_val, vecFl, vecCl, vecFp, vecCp, tmax, incr, b0, b0LF, delays, FieldList, nH, d1)
    
    
#Waiting time
    BeforeUp = L.make_delay(ThermalHF, bUp)
    AfterUp = L.make_delay(ThermalLF, aUp)
    
    BeforeDown = L.make_delay(ThermalLF, bDown)
    AfterDown = L.make_delay(ThermalHF, aDown)
    
#Shuttling
    #Selector
    SelectionShuttlingUp = np.zeros((16, 16), dtype=np.complex128)
    SelectionShuttlingUp[0, 0] = 1.0
    SelectionShuttlingUp[3, 3] = 1.0
    
    SelectionShuttlingDown = np.zeros((16, 16), dtype=np.complex128)
    SelectionShuttlingDown[0, 0] = 1.0
    SelectionShuttlingDown[3, 3] = 1.0
    SelectionShuttlingDown[6, 6] = 1.0
    SelectionShuttlingDown[15, 15] = 1.0
    
    
    
    Prop_Up = [[] for F in FieldList]
    Prop_Down = [[] for F in FieldList]
    for F in range(len(FieldList)):
        NewRates = R.compute_rates(FieldList[F], s2, tauc, tauf, CSA_val, vecFl, vecCl, vecFp = [1, 0, 0], vecCp = [1, 0, 0], lambda_f_ext=0.0, rho_f_ext=0.0, rH = nH/2.0)
        NewL = L.compute_free_precess_liouv(ChemShift_Val[0], ChemShift_Val[1], FieldList[F], NewRates, JCF)
        NewL_up =  L.ThermalCorr(NewL, NewRates)
        Prop_Up[F] = L.make_delay(NewL_up, incrShuttle)
        Prop_Down[F] = L.make_delay(NewL, incrShuttle)
        
    P_up = Prop_Up[0]
    P_down = Prop_Down[0]
    for F in range(len(FieldList)):
        P_up = Prop_Up[F] @ P_up
        P_down = P_down @ Prop_Down[F]
        
        
        
    
#t1 LF evolution
    P_t1 = []

    ta = T/2.0 - tau_180_C_LF/np.pi
    tb = 0.0 - tau_180_F_LF/np.pi - tau_180_C_LF/np.pi
    tc = T/2.0 - tau_180_F_LF/np.pi
    
    
    incr_a = incr1/2.
    incr_b = incr1/2. - T/(2.*(len(t1Delays)-1.))
    incr_c = T/(2.*(len(t1Delays)-1.))
    
    
    for t in range(len(t1Delays)):
        
        La = L.make_delay(ThermalLF, ta)
        Lb = L.make_delay(ThermalLF, tb)
        Lc = L.make_delay(ThermalLF, tc)
        
        P_t1.append(Lc @ P_180f_px_LF @ Lb @ P_180c_px_LF @ La)
        
        ta += incr_a
        tb += incr_b
        tc -= incr_c

        
    P_t1 = np.array(P_t1).reshape(-1, 1, 16, 16)
        
        
#t2 HF evolution
    P_t2 = []
    for t in t2Delays:
        P_t2.append(L.make_delay(ThermalHF, t))
    P_t2 = np.array(P_t2).reshape(-1, 16, 16)
    
    
    
    
######### PHASE CYCLING
    Phi1_x = [P_90f_px_LF, P_90f_mx_LF, P_90f_px_LF, P_90f_mx_LF, P_90f_px_LF, P_90f_mx_LF, P_90f_px_LF, P_90f_mx_LF]
    Phi1_y = [P_90f_py_LF, P_90f_my_LF, P_90f_py_LF, P_90f_my_LF, P_90f_py_LF, P_90f_my_LF, P_90f_py_LF, P_90f_my_LF]
    
    Phi2 = [P_90f_py_LF, P_90f_py_LF, P_90f_py_LF, P_90f_py_LF, P_90f_my_LF, P_90f_my_LF, P_90f_my_LF, P_90f_my_LF]
    
    
    Phi3 = [P_90c_mxpy, P_90c_mxpy, P_90c_mxmy, P_90c_mxmy, P_90c_mxpy, P_90c_mxpy, P_90c_mxmy, P_90c_mxmy]
    
    Receiver = [np.zeros((1, 16), dtype=np.complex128) for i in range(8)]
    Receiver[0][0, [4, 5]] = -1.0, -1j      #-x
    Receiver[1][0, [4, 5]] = 1.0, 1j      #x
    Receiver[2][0, [4, 5]] = 1.0, 1j      #x
    Receiver[3][0, [4, 5]] = -1.0, -1j      #-x
    Receiver[4][0, [4, 5]] = 1.0, 1j      #x
    Receiver[5][0, [4, 5]] = -1.0, -1j      #-x
    Receiver[6][0, [4, 5]] = -1.0, -1j      #-x
    Receiver[7][0, [4, 5]] = 1.0, 1j      #x
    
    
######### SIMULATE SEQUENCE
    M_LF = SelectionShuttlingUp @ AfterUp @ P_up @ BeforeUp @ P_steady
    
    
    
#Low field evolution
    M_PostT1_x = [Phi2[P] @ P_t1 @ Phi1_x[P] @ M_LF for P in range(len(Phi1_x))]
    M_PostT1_y = [Phi2[P] @ P_t1 @ Phi1_y[P] @ M_LF for P in range(len(Phi1_y))]
    
    
    PDOWN = P_180f_px @ AfterDown @ P_down @ SelectionShuttlingDown @ BeforeDown

    
#High field evolution
    Mdetect = [Receiver[S3E] @ P_t2 @ P_90c_py @ Phi3[S3E] @ P_S3E_HF @ P_180fc_px @ P_S3E_HF @ P_90c_px @ PDOWN for S3E in range(len(Phi3))]
    
    
    MDetect_x = [[] for i in range(len(Phi1_x))]
    MDetect_y = [[] for i in range(len(Phi1_x))]

    for P1 in range(len(Phi1_x)):
            MDetect_x[P1] = Mdetect[P1] @ M_PostT1_x[P1]
            MDetect_x[P1] = MDetect_x[P1][:, :, 0,0]
            
            MDetect_y[P1] = Mdetect[P1] @ M_PostT1_y[P1]
            MDetect_y[P1] = MDetect_y[P1][:, :, 0,0]
    
    
    return MDetect_x, MDetect_y
    
    
    
