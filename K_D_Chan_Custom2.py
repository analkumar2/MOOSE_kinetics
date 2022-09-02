# exec(open('../../Compilations/Kinetics/K_D_Chan_Custom2.py').read())
#KD channel taken from mod files of Migliore2018: kdb.mod
# Problems:

import numpy as np
import pickle
import pandas as pd
import moose
import matplotlib.pyplot as plt

#######################################################################
n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F = -0.033,-0.0087666, -3.28572168e-02, 3.39951757e-02, 0, 0, 1.65077372e-02, 4.08034957e-01

#######################################################################

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
    K_D.Ypower = 0.0
    K_D.Zpower = 0

    [nInf,nTau] = ChanGate(v,*[n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])

    xgate = moose.element( K_D.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = nInf/nTau
    xgate.tableB = 1.0/nTau

    return K_D


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