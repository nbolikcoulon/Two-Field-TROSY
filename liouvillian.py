# Imports
import numpy as np
from scipy import linalg

import sys


GAMMA_F = 2.51815e8
GAMMA_C = 6.72828e7



def ThermalCorr(liouvillian, rates):
    
    liouvillian_ = liouvillian.copy()
    
    liouvillian_[3, 0] -= 2.0 * rates['theta_f']
    liouvillian_[6, 0] -= 2.0 * rates['theta_c']
    liouvillian_[15, 0] -= 2.0 * rates['theta_cf']
    
    
    return liouvillian_


def compute_free_precess_liouv(wf, wc, b0, rates, J_CF):
    """Computes the relaxation matrix.

    The matrix is written in 16x16 basis, that is:
    E/2,
    {
    Fx, Fy, Fz,
    Cx, Cy, Cz,
    2FxCz, 2FyCz,
    2FzCx, 2FzCy,
    2FxCx, 2FxCy, 2FyCx, 2FyCy,
    2FzCz
    }

    Returns: numpy.matrix
    """
    

    liouvillian = np.zeros((16, 16))
    

    w_f = b0 * GAMMA_F * wf
    w_c = b0 * GAMMA_C * wc

    #
    # Populations - in units of gammaI*B0*hbar/(k*T)
#    liouvillian[3, 0] -= 2.0 * rates['theta_f']
#    liouvillian[6, 0] -= 2.0 * rates['theta_c']
#    liouvillian[15, 0] -= 2.0 * rates['theta_cf']
    liouvillian[3, 0] -= 0.0
    liouvillian[6, 0] -= 0.0
    liouvillian[15, 0] -= 0.0

    
    #
    # Chemical shift (I spin)
    # {Ix,Iy}
    liouvillian[1, 2] += w_f
    liouvillian[2, 1] -= w_f
    # {2IxSz,2IySz}
    liouvillian[7, 8] += w_f
    liouvillian[8, 7] -= w_f
    # {2IxSx,2IySx}
    liouvillian[11, 13] += w_f
    liouvillian[13, 11] -= w_f
    # {2IxSy,2IySy}
    liouvillian[12, 14] += w_f
    liouvillian[14, 12] -= w_f
    #
    # Chemical shift (S spin)
    # {Sx,Sy}
    liouvillian[4, 5] += w_c
    liouvillian[5, 4] -= w_c
    # {2IzSx,2IzSy}
    liouvillian[9, 10] += w_c
    liouvillian[10, 9] -= w_c
    # {2IxSx,2IxSy}
    liouvillian[11, 12] += w_c
    liouvillian[12, 11] -= w_c
    # {2IySx,2IySy}
    liouvillian[13, 14] += w_c
    liouvillian[14, 13] -= w_c
    #
    # Scalar coupling
    liouvillian[1, 8] += np.pi * J_CF
    liouvillian[8, 1] -= np.pi * J_CF
    liouvillian[2, 7] -= np.pi * J_CF
    liouvillian[7, 2] += np.pi * J_CF
    #
    liouvillian[4, 10] += np.pi * J_CF
    liouvillian[10, 4] -= np.pi * J_CF
    liouvillian[5, 9] -= np.pi * J_CF
    liouvillian[9, 5] += np.pi * J_CF
    #
    # Auto Relaxations
    liouvillian[1, 1] += rates['R2f']
    liouvillian[2, 2] += rates['R2f']
    liouvillian[3, 3] += rates['rho_f']
    liouvillian[4, 4] += rates['R2c']
    liouvillian[5, 5] += rates['R2c']
    liouvillian[6, 6] += rates['rho_c']
    liouvillian[7, 7] += rates['rhoa_f']
    liouvillian[8, 8] += rates['rhoa_f']
    liouvillian[9, 9] += rates['rhoa_c']
    liouvillian[10, 10] += rates['rhoa_c']
    liouvillian[11, 11] += rates['lambda_mq']
    liouvillian[12, 12] += rates['lambda_mq']
    liouvillian[13, 13] += rates['lambda_mq']
    liouvillian[14, 14] += rates['lambda_mq']
    liouvillian[15, 15] += rates['rho_2sp_cf']
    #
    # Longitudinal dipole-dipole cross-relaxation
    liouvillian[3, 6] += rates['sigma']
    liouvillian[6, 3] += rates['sigma']
    #
    # Longitudinal dipole-CSA cross-relaxation
    liouvillian[3, 15] += rates['delta_f']
    liouvillian[15, 3] += rates['delta_f']
    #
    liouvillian[6, 15] += rates['delta_c']
    liouvillian[15, 6] += rates['delta_c']
    #
    # Transverse dipole-CSA cross-relaxation
    liouvillian[1, 7] += rates['eta_f']
    liouvillian[7, 1] += rates['eta_f']
    liouvillian[2, 8] += rates['eta_f']
    liouvillian[8, 2] += rates['eta_f']
    #
    liouvillian[4, 9] += rates['eta_c']
    liouvillian[9, 4] += rates['eta_c']
    liouvillian[5, 10] += rates['eta_c']
    liouvillian[10, 5] += rates['eta_c']
    #
    # Multi-Quantum dipole-dipole cross-relaxation
    liouvillian[11, 14] -= rates['mu_mq']
    liouvillian[14, 11] -= rates['mu_mq']
    liouvillian[12, 13] += rates['mu_mq']
    liouvillian[13, 12] += rates['mu_mq']

    return liouvillian










def make_delay(liouvillian, time):
    return linalg.expm(-liouvillian * time)


def make_pulse_cf(liouvillian, w1_f=0.0, phase_f=0.0, w1_c=0.0, phase_c=0.0, pw=0.0):
    phase_rad_f = 0.5 * phase_f * np.pi
    phase_rad_c = 0.5 * phase_c * np.pi
    
    
    w1x_f = GAMMA_F * w1_f * np.cos(phase_rad_f)
    w1y_f = GAMMA_F * w1_f * np.sin(phase_rad_f)
    
    w1x_c = GAMMA_C * w1_c * np.cos(phase_rad_c)
    w1y_c = GAMMA_C * w1_c * np.sin(phase_rad_c)
    
    

    liouvillian_ = liouvillian.copy()
    
    # Pulses
    #

    liouvillian_[1, 3] -= w1y_f
    liouvillian_[3, 1] += w1y_f
    #
    liouvillian_[2, 3] += w1x_f
    liouvillian_[3, 2] -= w1x_f
    #
    liouvillian_[4, 6] -= w1y_c
    liouvillian_[6, 4] += w1y_c
    #
    liouvillian_[5, 6] += w1x_c
    liouvillian_[6, 5] -= w1x_c
    #
    liouvillian_[7, 11] += w1y_c
    liouvillian_[11, 7] -= w1y_c
    #
    liouvillian_[7, 12] -= w1x_c
    liouvillian_[12, 7] += w1x_c
    #
    liouvillian_[7, 15] -= w1y_f
    liouvillian_[15, 7] += w1y_f
    #
    liouvillian_[8, 13] += w1y_c
    liouvillian_[13, 8] -= w1y_c
    #
    liouvillian_[8, 14] -= w1x_c
    liouvillian_[14, 8] += w1x_c
    #
    liouvillian_[8, 15] += w1x_f
    liouvillian_[15, 8] -= w1x_f
    #
    liouvillian_[9, 11] += w1y_f
    liouvillian_[11, 9] -= w1y_f
    #
    liouvillian_[9, 13] -= w1x_f
    liouvillian_[13, 9] += w1x_f
    #
    liouvillian_[9, 15] -= w1y_c
    liouvillian_[15, 9] += w1y_c
    #
    liouvillian_[10, 12] += w1y_f
    liouvillian_[12, 10] -= w1y_f
    #
    liouvillian_[10, 14] -= w1x_f
    liouvillian_[14, 10] += w1x_f
    #
    liouvillian_[10, 15] += w1x_c
    liouvillian_[15, 10] -= w1x_c
    

    return linalg.expm(-liouvillian_ * pw)


def make_pulse_c(L, w1_c=0.0, phase_c=0.0, pw=0.0):
    return make_pulse_cf(L, w1_f=0.0, phase_f=0.0, w1_c=w1_c, phase_c=phase_c, pw=pw)

def make_pulse_f(L, w1_f=0.0, phase_f=0.0, pw=0.0):
    return make_pulse_cf(L, w1_f=w1_f, phase_f=phase_f, w1_c=0.0, phase_c=0.0, pw=pw)


def make_simpulse_cf(liouvillian, w1_f=0.0, phase_f=0.0, pw_f=0.0, w1_c=0.0, phase_c=0.0, pw_c=0.0):
    if pw_f < pw_c:

        pulse1 = make_pulse_c(liouvillian, w1_c=w1_c, phase_c=phase_c, pw=0.5 * (pw_c - pw_f))
        pulse2 = make_pulse_cf(liouvillian, w1_f=w1_f, phase_f=phase_f, w1_c=w1_c, phase_c=phase_c, pw=pw_f)

    else:

        pulse1 = make_pulse_f(liouvillian, w1_f=w1_f, phase_f=phase_f, pw=0.5 * (pw_f - pw_c))
        pulse2 = make_pulse_cf(liouvillian, w1_f=w1_f, phase_f=phase_f, w1_c=w1_c, phase_c=phase_c, pw=pw_c)

    return pulse1 @ pulse2 @ pulse1


def make_cf_pulse_from_file(liouvillian, w1_f=1.0, w1_c=1.0, pw=1.0, filename=None):
    pulse_phase_amplitudes = np.loadtxt(filename)
    pulse_phase_amplitudes[:, 1] /= 90.0
    pulse_phase_amplitudes[:, 3] /= 90.0
    pulse_phase_amplitudes[:, 0] *= w1_f
    pulse_phase_amplitudes[:, 2] *= w1_f

    propagator = np.eye(liouvillian.shape[0])
    pw_ = pw / float(len(pulse_phase_amplitudes))

    for ampl_f, phase_f, ampl_c, phase_c in pulse_phase_amplitudes:
        propagator = propagator.dot(
            make_pulse_cf(liouvillian, w1_f=ampl_f, phase_f=phase_f, w1_c=ampl_c, phase_c=phase_c, pw=pw_)
        )

    return propagator





def set_fzcz_eq():
    mag = np.zeros((16, 1))

    mag[0, 0] = 0.5
    mag[3, 0] = 1.0
    mag[6, 0] = 1.0

    return mag



def set_cz_eq():
    mag = np.zeros((16, 1))

    mag[0, 0] = 0.5
    mag[6, 0] = 1.0

    return mag


def set_cz_trosy_eq():
    mag = np.zeros((16, 1))

    mag[0, 0] = 0.5
    mag[6, 0] = 1.0
    mag[15, 0] = -1.0

    return 0.5 * mag


def set_2fzcz_eq():
    mag = np.zeros((16, 1))

    mag[0, 0] = 0.5
    mag[15, 0] = 1.0

    return mag


def set_fz_eq():
    mag = np.zeros((16, 1))

    mag[0, 0] = 0.5
    mag[3, 0] = 1.0

    return mag


def get_cz(mag):
    mag_cz = mag[6, 0]

    return mag_cz


def get_cz_trosy(mag):
    mag_czTrosy = mag[6, 0] - mag[15, 0]

    return mag_czTrosy


def get_fz(mag):
    mag_fz = mag[1, 0]

    return mag_fz
