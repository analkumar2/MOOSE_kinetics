# exec(open('../../Compilations/Kinetics/K_31_Chan_Hay2011.py').read())
# Base is Nandi 2022. They took it from Hay 2011. Not Parameterized. 

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

def K_31_Chan(name):
    K_31 = moose.HHChannel( '/library/' + name )
    K_31.Ek = EK
    K_31.Gbar = 300.0*SOMA_A
    K_31.Gk = 0.0
    K_31.Xpower = 1.0
    K_31.Ypower = 0
    K_31.Zpower = 0

    vshift = 0

    mInf =  1/(1+np.exp(((v -(18.700 + vshift))/(-9.700))))
    mTau =  0.2*20.000/(1+np.exp(((v -(-46.560 + vshift))/(-44.140))))

    xgate = moose.element( K_31.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = mInf/mTau
    xgate.tableB = 1.0/mTau

    return K_31


if __name__ == "__main__":
    moose.Neutral('library')
    K_31_Chan('K_31_Chan')
    plt.figure()
    plt.plot(v, moose.element('library/K_31_Chan/gateX').tableA/moose.element('library/K_31_Chan/gateX').tableB, label='nInf')    
    plt.ylabel('Inf')
    plt.legend()
    plt.figure()
    plt.plot(v, 1/moose.element('library/K_31_Chan/gateX').tableB, label='nTau')
    plt.ylabel('Tau')
    plt.legend()
    plt.show()