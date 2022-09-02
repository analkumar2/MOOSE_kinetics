# K_Tst channel K_Tst.mod Hay2011
# exec(open('K_Tst_Chan_(Hay2011).py').read())

import numpy as np
import pickle
import pandas as pd
import moose

SOMA_A = 3.14e-8
F = 96485.3329
R = 8.314
celsius = 32
dt = 0.05e-3
ENa = 0.050
EK = -0.085
Eh = -0.045
ECa = 0.128
Em = -0.090

Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
v = np.linspace(Vmin,Vmax, Vdivs)
Camin = 0.01e-3
Camax = 1e-3
Cadivs = 4000
ca = np.linspace(Camin,Camax, Cadivs)

def K_Tst_Chan(name):
    K_Tst = moose.HHChannel( '/library/' + name )
    K_Tst.Ek = EK
    K_Tst.Gbar = 300.0*SOMA_A
    K_Tst.Gk = 0.0
    K_Tst.Xpower = 4.0
    K_Tst.Ypower = 1.0
    K_Tst.Zpower = 0.0

    qt = 2.3**((34-21)/10)
    V = v + 0.010
    mInf =  1/(1 + np.exp(-(V*1e3+0)/19))
    mTau =  (0.34+0.92*np.exp(-((V*1e3+71)/59)**2))/qt
    hInf =  1/(1 + np.exp(-(V*1e3+66)/-10))
    hTau =  (8+49*np.exp(-((V*1e3+73)/23)**2))/qt

    xgate = moose.element( K_Tst.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = mInf/mTau*1e3
    xgate.tableB = 1/mTau*1e3

    ygate = moose.element( K_Tst.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hInf/hTau*1e3
    ygate.tableB = 1/hTau*1e3

    return K_Tst
