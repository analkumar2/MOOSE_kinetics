# exec(open('../../Compilations/Kinetics/K_custom1_Chan_Custom3.py').read())
# K_custom1 channel parameterized kinetics.

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
l_vhalf_inf, l_slope_inf, l_A, l_B, l_C, l_D, l_E, l_F = 0.013,-0.0087666, 1.26E-02,1.73E-02,0.00E+00,0.00E+00,3.43E-02,1.02E-01
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

def K_custom1_Chan(name):
    K_custom1 = moose.HHChannel( '/library/' + name )
    K_custom1.Ek = EK
    K_custom1.Gbar = 300.0*SOMA_A
    K_custom1.Gk = 0.0
    K_custom1.Xpower = 1.0
    K_custom1.Ypower = 1.0
    K_custom1.Zpower = 0

    [nInf,nTau] = ChanGate(v,*[n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])
    [lInf,lTau] = ChanGate(v,*[l_vhalf_inf, l_slope_inf, l_A, l_B, l_C, l_D, l_E, l_F])

    xgate = moose.element( K_custom1.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = nInf/nTau
    xgate.tableB = 1.0/nTau

    ygate = moose.element( K_custom1.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = lInf/lTau
    ygate.tableB = 1.0/lTau

    return K_custom1


if __name__ == "__main__":
    [nInf,nTau] = ChanGate(v,*[n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])
    [lInf,lTau] = ChanGate(v,*[l_vhalf_inf, l_slope_inf, l_A, l_B, l_C, l_D, l_E, l_F])
    plt.figure()
    plt.plot(v, nInf, label='nInf')
    plt.plot(v, lInf, label='lInf')
    plt.ylabel('Inf')
    plt.legend()
    plt.figure()
    plt.plot(v, nTau, label='nTau')
    plt.plot(v, lTau, label='lTau')
    plt.ylabel('Tau')
    plt.legend()
    plt.show()