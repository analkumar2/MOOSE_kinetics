# exec(open('../../Compilations/Kinetics/K_DR_Chan_Custom3.py').read())
# K_DR channel parameterized kinetics. Base is experimental kinetics by Migliore2018. 4 gates instead of 1

import numpy as np
import pickle
import pandas as pd
import moose
import matplotlib.pyplot as plt

SOMA_A = 3.14e-8
F = 96485.3329
R = 8.314
celsius = 32
dt = 0.05e-3
ENa = 0.092
EK = -0.099
Eh = -0.030
ECa = 0.140
Em = -0.065

#################################
n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F = 0.013,0.0087666, 1.26E-02,1.73E-02,0.00E+00,0.00E+00,3.43E-02,1.02E-01
#################################

Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
# dV = (Vmax-Vmin)/Vdivs
# v = np.arange(Vmin,Vmax, dV)
v = np.linspace(Vmin,Vmax, Vdivs)
Camin = 0.04e-3
Camax = 1
Cadivs = 8000
# dCa = (Camax-Camin)/Cadivs
# ca = np.arange(Camin,Camax, dCa)
ca = np.linspace(Camin,Camax, Cadivs)

def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
    # alge model
    Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
    yl = (v-A)/-B
    yr = (v-A)/E
    Tau = (C + (1 + yl/(np.sqrt(1+yl**2)))/2) * (D + (1 + yr/(np.sqrt(1+yr**2)))/2) * F
    Tau[Tau<0.00002] = 0.00002
    return [Inf,Tau]

def K_DR_Chan(name):
    K_DR = moose.HHChannel( '/library/' + name )
    K_DR.Ek = EK
    K_DR.Gbar = 300.0*SOMA_A
    K_DR.Gk = 0.0
    K_DR.Xpower = 4.0
    K_DR.Ypower = 0.0
    K_DR.Zpower = 0

    [nInf,nTau] = ChanGate(v,*[n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])

    xgate = moose.element( K_DR.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = nInf/nTau
    xgate.tableB = 1.0/nTau

    return K_DR


if __name__ == "__main__":
    [nInf,nTau] = ChanGate(v,*[n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])
    plt.figure()
    plt.plot(v, nInf, label='nInf')
    plt.ylabel('Inf')
    plt.legend()
    plt.figure()
    plt.plot(v, nTau, label='nTau')
    plt.ylabel('Tau')
    plt.legend()
    plt.show()