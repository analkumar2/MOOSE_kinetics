# exec(open('../../Compilations/Kinetics/Na_T_Chan_Hay2011.py').read())
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
ENa = 53e-3 #0.092 #from Deepanjali data
EK = -107e-3 #-0.099 #from Deepanjali data
Eh = -0.030
ECa = 0.140 #from Deepanjali data
Em = -0.065


#################################
n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F = -0.04230235,0.006, -0.04476696,  0.02195918,0.01651612,  0.04558796,0.01921329,  0.00093119
l_vhalf_inf, l_slope_inf, l_A, l_B, l_C, l_D, l_E, l_F = -0.066,-0.006, -0.07071569,  0.02638251,0.01824087,  0.06117709,0.01613637,  0.00866347
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

def Na_T_Chan(name):
    Na_T = moose.HHChannel( '/library/' + name )
    Na_T.Ek = ENa
    Na_T.Gbar = 300.0*SOMA_A
    Na_T.Gk = 0.0
    Na_T.Xpower = 3.0
    Na_T.Ypower = 1.0
    Na_T.Zpower = 0

    [nInf,nTau] = ChanGate(v,*[n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])
    [lInf,lTau] = ChanGate(v,*[l_vhalf_inf, l_slope_inf, l_A, l_B, l_C, l_D, l_E, l_F])
    


    xgate = moose.element( Na_T.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = nInf/nTau
    xgate.tableB = 1.0/nTau

    ygate = moose.element( Na_T.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = lInf/lTau
    ygate.tableB = 1.0/lTau

    return Na_T


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