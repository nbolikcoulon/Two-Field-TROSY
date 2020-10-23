#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 18:30:47 2019

@author: nbc
"""

import pulsesequence as PS
import fitfunctions as FitF
import ratetrosy as R
#import rate as R


import numpy as np

from matplotlib import pyplot as plt


import nmrglue as ng


import sys



Navy = (0.0, 0.0, 0.501961)
Deeppink = (1.0, 0.0784314, 0.576471)
DarkCyan = (0.0, 0.545098, 0.545098)
Cornflowerblue = (0.392157, 0.584314, 0.929412)
BlueViolet = (0.541176, 0.168627, 0.886275)
Mediumvioletred = (0.780392, 0.0823529, 0.521569)
DarkOrange = (1.0, 0.54902, 0.0)
Crimson = (0.862745, 0.0784314, 0.235294)
Mediumseegreen = (0.235294, 0.701961, 0.443137)
Indigo = (0.294118, 0.0, 0.509804)
Darkgreen = (0.0, 0.392157, 0.0)
blue = (0.0, 0.0, 1.0)
        
colors = [Indigo, Mediumseegreen, Deeppink, Cornflowerblue, Mediumvioletred, DarkOrange, Crimson, Darkgreen, BlueViolet, Mediumvioletred, blue, Navy, DarkCyan]
for i in np.linspace(0, 1, 100):
    colors.append((i, i, i))



GAMMA_F = 2.51815e8
GAMMA_C = 6.72828e7



s2 = 0.8
tauf = 100e-12


t1ScaleRange = np.arange(0.5, 3.1, 0.25)



#tc = 25 ns
tauc = 25e-9
b0 = [21.15, 14.1]
b0LF = [2.5, 3.0]
d1 = [0.53, 1.11]


##tc = 100 ns
#tauc = 100e-9
#b0 = [21.15, 14.1]
#b0LF = [2.5, 3.0]
#d1 = [0.78, 1.90]




Molecules = ["3FTyr", "4FPhe"]


CSA_F_XX = [293e-6, 300e-6]
CSA_F_YY = [250e-6, 219e-6]
CSA_F_ZZ = [407e-6, 358e-6]

CSA_C_XX = [-47e-6, -82e-6]
CSA_C_YY = [28e-6, 29e-6]
CSA_C_ZZ = [82e-6, 80e-6]


Angle_F_XX_ZZ = [np.deg2rad(17.0), np.deg2rad(0.0)]
Angle_C_XX_ZZ = [np.deg2rad(16.0), np.deg2rad(0.0)]


CSA_val = []
vecFl = []
vecCl = []
ChemShift_val = np.zeros((len(CSA_F_XX), 2))
for i in range(len(CSA_F_XX)):
    CSA_F_long = CSA_F_XX[i] - CSA_F_YY[i]
    CSA_F_perp = CSA_F_ZZ[i] - CSA_F_YY[i]
    CSA_C_long = CSA_C_XX[i] - CSA_C_YY[i]
    CSA_C_perp = CSA_C_ZZ[i] - CSA_C_YY[i]
    CSA_val.append([CSA_F_long, CSA_F_perp, CSA_C_long, CSA_C_perp])
    
    ChemShift_F = np.average([CSA_F_XX[i], CSA_F_YY[i], CSA_F_ZZ[i]])
    ChemShift_C = np.average([CSA_C_XX[i], CSA_C_YY[i], CSA_C_ZZ[i]])
    ChemShift_val[i, 0] = ChemShift_F
    ChemShift_val[i, 1] = ChemShift_C

    vecFl.append([0.0, -np.sin(Angle_F_XX_ZZ[i]), np.cos(Angle_F_XX_ZZ[i])])
    vecCl.append([0.0, -np.sin(Angle_C_XX_ZZ[i]), np.cos(Angle_C_XX_ZZ[i])])
    

    
nH = [0.0, 0.0]
    

delays = [25e-3, 40e-3, 5e-3, 350e-3, 1e-3] #bUp, aUp, bDown, aDown, incr
ShuttlingTime = 100e-3
    
    
#Fit field profile
FieldCalibration = open("FieldCalibration.txt", "r")
#This is a file containing the height and the corresponding field

with FieldCalibration as input:
    FC = list(zip(*(line.strip().split("\t") for line in input)))
FC = [[float(FC[col][line]) for col in range(2)] for line in range(len(FC[0]))]   #it is a one column vector containing the height and the field at each position
heights = [FC[col][0] for col in range(len(FC))]
fields = [FC[col][1] for col in range(len(FC))]



xlim = [[[321, 312], [26, 16]], [[297, 288], [14, 4]]]




def CalcNPoints(B0HF, B0LF, tmax):
    SWppm_F = 10e-6
    SWppm_C = 10e-6
    N_F = 1. + SWppm_F*tmax[0]*GAMMA_F*B0LF/(2.*np.pi)
    N_C = 1. + SWppm_C*tmax[1]*GAMMA_C*B0HF/(2.*np.pi)
    
    return [int(N_F), int(N_C)]



                
b0Ref = 14.1
                
FileName = 'D1_1FST2PT_3FTyr_tc25.0ns.txt'
File = open(FileName, 'r')
File.readline()
check = 0
while check == 0:
    Line = File.readline()
    Line = Line.split('\n')
    Line = Line[0]
    Line = Line.split('\t')
    
    if round(float(Line[0]), 2) == b0Ref:
        d1Ref = float(Line[2])
        check = 1
        
RatesHF = R.compute_rates(b0Ref, s2, 25e-9, tauf, CSA_val[0], vecFl[0], vecCl[0], vecFp = [1, 0, 0], vecCp = [1, 0, 0], lambda_f_ext=0.0, rho_f_ext=0.0, rH = 0.0, Ftrosy = 1.0)
t2MaxRef = 3.0/RatesHF['TrosyC']
t1MaxRef = 1.25/RatesHF['TrosyF']

NPRef = CalcNPoints(b0Ref, b0Ref, [t1MaxRef, t2MaxRef])
incrRef = [t1MaxRef/(NPRef[0]-1.), t2MaxRef/(NPRef[1]-1.)]
t1DelaysRef = np.arange(0.0, t1MaxRef+incrRef[0], incrRef[0])
        
tRef = 0.0
tPS = 1./240.
time = tPS + t2MaxRef + d1Ref
for n in range(NPRef[0]):
    tRef += time + t1DelaysRef[n]
tRef = tRef*8.
        

        
#################################




File = open('ExpTime_tc' + str(int(1e9*tauc))  + 'ns.txt',  'w')
File.write('Molecule\ttauc\tt1 scale\tExp  time (s)\tExp  time (min)')


for Mol in range(len(Molecules)):
    
    scalingField = b0[Mol]/fields[0]
    ScaleField = [scalingField*F for F in fields]
    HigherCoefs, MiddleCoefs, LowerCoefs = FitF.FitField(heights, ScaleField)
    ShuttledHeight = heights[ScaleField.index(min(ScaleField, key=lambda x:abs(x-b0LF[Mol])))]
    FieldList = FitF.FieldList(ShuttlingTime, ShuttledHeight, delays[-1], LowerCoefs, MiddleCoefs, HigherCoefs)
    
    RHF =  R.compute_rates(b0[Mol], s2, tauc, tauf, CSA_val[Mol], vecFl[Mol], vecCl[Mol], vecFp = [1, 0, 0], vecCp = [1, 0, 0], lambda_f_ext=0.0, rho_f_ext=0.0, rH = nH[Mol]/2.0)
    RLF =  R.compute_rates(b0LF[Mol], s2, tauc, tauf, CSA_val[Mol], vecFl[Mol], vecCl[Mol], vecFp = [1, 0, 0], vecCp = [1, 0, 0], lambda_f_ext=0.0, rho_f_ext=0.0, rH = nH[Mol]/2.0)
    
    AvT2F_LF = 1./RLF['R2f']
    AvT2TrosyC_HF = 1./RHF['TrosyC']
    
    File.write('\n')
    
    for t1Scale in range(len(t1ScaleRange)):
        
        tmax = [t1ScaleRange[t1Scale] * AvT2F_LF, 3.0*AvT2TrosyC_HF]

        NP = CalcNPoints(b0[Mol], b0LF[Mol], tmax)
        incr = [tmax[0]/(NP[0]-1.), tmax[1]/(NP[1]-1.)]
        
        
        t1Delays = np.arange(0.0, tmax[0]+incr[0], incr[0])
        exptime = 0.0
        tPS = sum(delays) - 1e-3 + 2.*ShuttlingTime + 3./(4.*240.)
        time = tPS + tmax[1] + d1[Mol]
        for n in range(NP[0]):
            exptime += time + t1Delays[n]
        exptime = exptime*16.
        print(Molecules[Mol], tauc, t1ScaleRange[t1Scale], exptime/60.)
        
        File.write('\n')
        File.write(Molecules[Mol])
        File.write('\t')
        File.write(str(tauc))
        File.write('\t')
        File.write(str(t1ScaleRange[t1Scale]))
        File.write('\t')
        File.write(str(exptime))
        File.write('\t')
        File.write(str(exptime/60.))

        
        
        ScaleInt = 80. * b0[Mol]/b0Ref * tRef/exptime
        
        
        SWHz_F = (NP[0]-1.)/tmax[0]
        SWHz_C = (NP[1]-1.)/tmax[1]
        
                    
        ChemShift_val_OFFSET = [0.0, 0.0]
            
            
        yT1, xT1 = PS.CdetecTROSY_S3E_ST2PT(s2, tauc, tauf, ChemShift_val_OFFSET, CSA_val[Mol], vecFl[Mol], vecCl[Mol], [1, 0, 0], [1, 0, 0], tmax, incr, b0[Mol], b0LF[Mol], delays, FieldList, nH[Mol], d1[Mol])
        
        #Reconstruct direct dimension
        yT1 = sum(yT1[i] for i in range(len(yT1)))
        xT1 = sum(xT1[i] for i in range(len(xT1)))
        
        
        SignalY = ScaleInt * yT1
        SignalX = ScaleInt * xT1
            
        NoiseEcho = [np.sqrt(ScaleInt) * np.random.normal(loc=0.0, scale=2, size=np.asarray(SignalY[i]).shape) for i in range(len(SignalY))]
        NoiselAntiEcho = [np.sqrt(ScaleInt) * np.random.normal(loc=0.0, scale=2, size=np.asarray(SignalX[i]).shape) for i in range(len(SignalX))]
            
        SignalY = [SignalY[i] + NoiseEcho[i] for i in range(len(SignalY))]
        SignalX = [SignalX[i] + NoiselAntiEcho[i] for i in range(len(SignalX))]
        SignalY_N = [NoiseEcho[i] for i in range(len(NoiseEcho))]
        SignalX_N = [NoiselAntiEcho[i] for i in range(len(NoiselAntiEcho))]
        
            
        
            
            
        #FT in direct dimension        
        SignalY = ng.process.proc_base.sp(np.asarray(SignalY), off=0.5, end=0.98, pow=2.0)
        SignalY = ng.process.proc_base.zf_double(SignalY, 1)
        SignalY_N = ng.process.proc_base.sp(np.asarray(SignalY_N), off=0.5, end=0.98, pow=2.0)
        SignalY_N = ng.process.proc_base.zf_double(SignalY_N, 1)
        
        SignalX = ng.process.proc_base.sp(np.asarray(SignalX), off=0.5, end=0.98, pow=2.0)
        SignalX = ng.process.proc_base.zf_double(SignalX, 1)
        SignalX_N = ng.process.proc_base.sp(np.asarray(SignalX_N), off=0.5, end=0.98, pow=2.0)
        SignalX_N = ng.process.proc_base.zf_double(SignalX_N, 1)
        
        SignalY[0, :] /= 2.0
        SignalX[0, :] /= 2.0
        SignalY_N[0, :] /= 2.0
        SignalX_N[0, :] /= 2.0
        
        
        SignalY_FT = ng.process.proc_base.fft(SignalY)
        SignalX_FT = ng.process.proc_base.fft(SignalX)
        SignalY_N_FT = ng.process.proc_base.fft(SignalY_N)
        SignalX_N_FT = ng.process.proc_base.fft(SignalX_N)
        
            
        #Phasing in direct dimension
        #    p0, p1 = ng.process.proc_autophase.manual_ps(SignalY_FT[0])
        #    sys.exit()
        p0, p1 = 0.0, 0.0
        SignalY_FTphased = ng.proc_base.ps(SignalY_FT, p0=p0, p1=p1)
        SignalX_FTphased = ng.proc_base.ps(SignalX_FT, p0=p0, p1=p1)
        SignalY_N_FTphased = ng.proc_base.ps(SignalY_N_FT, p0=p0, p1=p1)
        SignalX_N_FTphased = ng.proc_base.ps(SignalX_N_FT, p0=p0, p1=p1)
        
        
            
        #Transpose
        TransposeY = np.transpose(SignalY_FTphased)
        TransposeX = np.transpose(SignalX_FTphased)
        TransposeY_N = np.transpose(SignalY_N_FTphased)
        TransposeX_N = np.transpose(SignalX_N_FTphased)
            
            
        #Reconstruct indirect dimension
        Reconstruct = TransposeX.real + 1j*TransposeY.real
        Reconstruct_N = TransposeX_N.real + 1j*TransposeY_N.real
        
        
        #FT in indirect dimension
        Reconstruct = ng.process.proc_base.sp(Reconstruct, off=0.5, end=0.98, pow=2.0)
        Reconstruct = ng.process.proc_base.zf_double(Reconstruct, 1)
        Spectrum = ng.process.proc_base.fft(Reconstruct)
        
        Reconstruct_N = ng.process.proc_base.sp(Reconstruct_N, off=0.5, end=0.98, pow=2.0)
        Reconstruct_N = ng.process.proc_base.zf_double(Reconstruct_N, 1)
        Spectrum_N = ng.process.proc_base.fft(Reconstruct_N)
        
        #Phasing in indirect dimension
        #    p0, p1 = ng.process.proc_autophase.manual_ps(Spectrum[0])
        p0, p1 = 90.0, 0.0
        Spectrum_Phased = ng.proc_base.ps(Spectrum, p0=p0, p1=p1)
        Spectrum_N_Phased = ng.proc_base.ps(Spectrum_N, p0=p0, p1=p1)
        
        Spectrum_Final = np.transpose(Spectrum_Phased.real)
        Spectrum_N_Final = np.transpose(Spectrum_N_Phased.real)
            
            
                
        obsF = b0LF[Mol]*GAMMA_F*1e-6/(2.0*np.pi)
        obsC = b0[Mol]*GAMMA_C*1e-6/(2.0*np.pi)
            
        carF = ChemShift_val[Mol, 0] * GAMMA_F*b0LF[Mol]/(2.0*np.pi)
        carC = ChemShift_val[Mol, 1] * GAMMA_C*b0[Mol]/(2.0*np.pi)
            
        
        NP = Spectrum_Final.shape
        
            
            
        udic = {
            'ndim': 2,
            0: {'car': carF,
                'complex': False,
                'encoding': 'Rance-Kay',
                'freq': True,
                'label': '19F',  
                'obs': obsF,
                'size': NP[0],
                'sw': SWHz_F,
                'time': False},
            1: {'car': carC,
                'complex': False,
                'encoding': 'direct',
                'freq': True,
                'label': '13C',
                'obs': obsC,
                'size': NP[1],
                'sw': SWHz_C,
                'time': False}
                }
            
            
            
            
        dic = ng.sparky.create_dic(udic)
        
        name = "Spectra/" + str(int(1e9*tauc)) + "ns/" + Molecules[Mol] + "/2FS3E_" + Molecules[Mol] + "_" + str(round(b0LF[Mol], 2)) + "_" + str(round(b0[Mol], 2)) + "_tc" + str(tauc) + "_" + str(t1ScaleRange[t1Scale]) + "t1.ucsf"
        name_N = "Spectra/" + str(int(1e9*tauc)) + "ns/" + Molecules[Mol] + "/Noise/Noise2FS3E_" + Molecules[Mol] + "_" + str(round(b0LF[Mol], 2)) + "_" + str(round(b0[Mol], 2)) + "_tc" + str(tauc)  + "_" + str(t1ScaleRange[t1Scale]) + "t1.ucsf"
                
        ng.sparky.write(name, dic, np.asarray(Spectrum_Final).astype('float32'), overwrite=True)
        ng.sparky.write(name_N, dic, np.asarray(Spectrum_N_Final).astype('float32'), overwrite=True)
        
    
        dic_2F, data_2F = ng.sparky.read(name)
        udic_2F = ng.sparky.guess_udic(dic_2F, data_2F)
        
        uc_C_2F = ng.sparky.make_uc(dic_2F, data_2F, dim=1)
        ppm_C_2F = uc_C_2F.ppm_scale()
        uc_F_2F = ng.sparky.make_uc(dic_2F, data_2F, dim=0)
        ppm_F_2F = uc_F_2F.ppm_scale()
        
        Intensity_2F = max([max(data_2F[i]) for i in range(len(data_2F))])
        pos_Max = np.where(data_2F == Intensity_2F)
        pos_Max = [pos_Max[0][0], pos_Max[1][0]]
        
        
        FcrossSection_2F = data_2F[:, pos_Max[1]]
        CcrossSection_2F = data_2F[pos_Max[0], :]
        
        
        
        
        
        
        dic_Noise, data_Noise = ng.sparky.read(name_N)
        STDNoise = np.std(data_Noise)
        SN = Intensity_2F/STDNoise
        
        
#        print(STDNoise, SN)
        
        
        
        
        plt.plot(ppm_F_2F, FcrossSection_2F.real, label='S/N: ' + str(round(SN, 2)))
        plt.title(Molecules[Mol] + ' - tc = ' + str(tauc) + ' ns - ' + str(t1ScaleRange[t1Scale]) + ' x t1')
        plt.xlabel('19F - ppm')
        plt.xlim(xlim[Mol][0])
        plt.legend(loc='best')
        plt.savefig('Fig/' + str(int(1e9*tauc)) + "ns/" + Molecules[Mol] + '/19F/19F_2FS3E_' + str(t1ScaleRange[t1Scale]) + 't1.pdf', format='pdf')
        plt.close()
        
        
        plt.plot(ppm_C_2F, CcrossSection_2F.real, label='S/N: ' + str(round(SN, 2)))
        plt.title(Molecules[Mol] + ' - tc = ' + str(tauc) + ' ns - ' + str(t1ScaleRange[t1Scale]) + ' x t1')
        plt.xlabel('13C - ppm')
        plt.xlim(xlim[Mol][1])
        plt.legend(loc='best')
        plt.savefig('Fig/' + str(int(1e9*tauc)) + "ns/" + Molecules[Mol] + '/13C/13C_2FS3E_' + str(t1ScaleRange[t1Scale]) + 't1.pdf', format='pdf')
        plt.close()