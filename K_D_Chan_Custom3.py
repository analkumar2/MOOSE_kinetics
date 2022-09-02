# exec(open('../../Compilations/Kinetics/K_D_Chan_Custom3.py').read())
# Not based on any exp. Copy of K_A but inactivation at ~100ms tau

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
n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F = 0.0112,0.017, -8.78e-3,5.63e-2,0,0,2.65e-2,1.05e-2
l_vhalf_inf, l_slope_inf, l_min, l_m, l_cm = -0.056,-0.00877, 0.002, 0.26, 0.50
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

def K_D_Chan(name):
    K_D = moose.HHChannel( '/library/' + name )
    K_D.Ek = EK
    K_D.Gbar = 300.0*SOMA_A
    K_D.Gk = 0.0
    K_D.Xpower = 1.0
    K_D.Ypower = 1.0
    K_D.Zpower = 0

    [nInf,nTau] = ChanGate(v,*[n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])
    lInf = ChanGate(v,*[l_vhalf_inf, l_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])[0]


    lTau = l_m*(v+l_cm)
    lTau[lTau<l_min] = l_min


    xgate = moose.element( K_D.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = nInf/nTau
    xgate.tableB = 1.0/nTau

    ygate = moose.element( K_D.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = lInf/lTau
    ygate.tableB = 1.0/lTau

    return K_D


if __name__ == "__main__":
    [nInf,nTau] = ChanGate(v,*[n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])
    lInf = ChanGate(v,*[l_vhalf_inf, l_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])[0]
    lTau = l_m*(v+l_cm)
    lTau[lTau<l_min] = l_min

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