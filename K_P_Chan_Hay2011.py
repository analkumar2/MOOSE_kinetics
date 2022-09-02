# exec(open('../../Compilations/Kinetics/K_P_Chan_Hay2011.py').read())
# Base is Nandi 2022. They took it from Hay 2011. Parameterized. 

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
ENa = 0.092 #from Deepanjali data
EK = -0.099 #from Deepanjali data
Eh = -0.030
ECa = 0.140 #from Deepanjali data
Em = -0.065


#################################
n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F = -0.0143,0.0146, -5.79385816e-02,  4.05353979e-02,4.64575934e-05,  8.25998547e-02,1.96014286e-02,  5.71474459e-02
l_vhalf_inf, l_slope_inf, l_A, l_B, l_C, l_D, l_E, l_F = -0.054,-0.011, -5.29123735e-02,  3.45061811e-02,5.91067020e-02,  6.40390236e-10,3.36911564e-02,  1.79793456e+00
#################################

Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
# dV = (Vmax-Vmin)/Vdivs
# v = np.arange(Vmin,Vmax, dV)
v = np.linspace(Vmin,Vmax, Vdivs)
Camin = 1e-12
Camax = 3
Cadivs = 4000
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

def K_P_Chan(name):
    K_P = moose.HHChannel( '/library/' + name )
    K_P.Ek = EK
    K_P.Gbar = 300.0*SOMA_A
    K_P.Gk = 0.0
    K_P.Xpower = 2.0
    K_P.Ypower = 1.0
    K_P.Zpower = 0

    [nInf,nTau] = ChanGate(v,*[n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])
    [lInf,lTau] = ChanGate(v,*[l_vhalf_inf, l_slope_inf, l_A, l_B, l_C, l_D, l_E, l_F])
    


    xgate = moose.element( K_P.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = nInf/nTau
    xgate.tableB = 1.0/nTau

    ygate = moose.element( K_P.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = lInf/lTau
    ygate.tableB = 1.0/lTau

    return K_P


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