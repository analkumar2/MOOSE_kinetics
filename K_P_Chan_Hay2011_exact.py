# exec(open('../../Compilations/Kinetics/K_P_Chan_Hay2011.py').read())
# Base is Nandi 2022. They took it from Hay 2011. Not parameterized. 

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


Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
# dV = (Vmax-Vmin)/Vdivs
# v = np.arange(Vmin,Vmax, dV)
v = np.linspace(Vmin,Vmax, Vdivs)*1e3 #1e3 because this is converted from NEURON which uses mV and ms
Camin = 1e-12
Camax = 3
Cadivs = 4000
# dCa = (Camax-Camin)/Cadivs
# ca = np.arange(Camin,Camax, dCa)
ca = np.linspace(Camin,Camax, Cadivs)

def K_P_Chan(name):
    K_P = moose.HHChannel( '/library/' + name )
    K_P.Ek = EK
    K_P.Gbar = 300.0*SOMA_A
    K_P.Gk = 0.0
    K_P.Xpower = 2.0
    K_P.Ypower = 1.0
    K_P.Zpower = 0

    vshift = 0
    tauF = 1

    def vtrap(x,y):
        if abs(x/y)<1e-6:
            return y * (1 - x / y / 2)
        else:
            return x / (np.exp(x / y) - 1)

    qt = 2.3**((celsius-21)/10)
    mInf =  1 / (1 + np.exp(-(v - (-14.3 + vshift)) / 14.6))
    mTau = tauF * (1.25+13*np.exp(-(v - vshift) * 0.026))/qt
    mTau[v < -50 + vshift] = tauF * (1.25+175.03*np.exp(-(v[v < -50 + vshift] - vshift) * -0.026))/qt
    hInf =  1/(1 + np.exp(-(v - (-54 + vshift))/-11))
    hTau =  (360+(1010+24*(v - (-55 + vshift)))*np.exp(-((v - (-75 + vshift))/48)**2))/qt

    xgate = moose.element( K_P.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = mInf/mTau*1e3
    xgate.tableB = 1.0/mTau*1e3

    ygate = moose.element( K_P.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hInf/hTau*1e3
    ygate.tableB = 1.0/hTau*1e3

    return K_P


if __name__ == "__main__":
    moose.Neutral('library')
    K_P_Chan('K_P_Chan')
    plt.figure()
    plt.plot(v, moose.element('library/K_P_Chan/gateX').tableA/moose.element('library/K_P_Chan/gateX').tableB, label='nInf')
    plt.plot(v, moose.element('library/K_P_Chan/gateY').tableA/moose.element('library/K_P_Chan/gateY').tableB, label='lInf')
    plt.ylabel('Inf')
    plt.legend()
    plt.figure()
    plt.plot(v, 1/moose.element('library/K_P_Chan/gateX').tableB, label='nTau')
    plt.plot(v, 1/moose.element('library/K_P_Chan/gateY').tableB, label='lTau')
    plt.ylabel('Tau')
    plt.legend()
    plt.show()