# K_Pst channel K_Pst.mod Hay2011
# exec(open('K_Pst_Chan_(Hay2011).py').read())

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

def K_Pst_Chan(name):
    K_Pst = moose.HHChannel( '/library/' + name )
    K_Pst.Ek = EK
    K_Pst.Gbar = 300.0*SOMA_A
    K_Pst.Gk = 0.0
    K_Pst.Xpower = 2.0
    K_Pst.Ypower = 1.0
    K_Pst.Zpower = 0.0

    qt = 2.3**((34-21)/10)
    V = v + 0.010
    mInf =  (1/(1 + np.exp(-(V*1e3+1)/12)))
    VV = V[V<-0.050]
    mTau =  (1.25+175.03*np.exp(-VV*1e3 * -0.026))/qt
    VV = V[V>=-0.050]
    mTau = np.append(mTau, ((1.25+13*np.exp(-VV*1e3*0.026)))/qt)
    hInf =  1/(1 + np.exp(-(V*1e3+54)/-11))
    hTau =  (360+(1010+24*(V*1e3+55))*np.exp(-((V*1e3+75)/48)**2))/qt

    xgate = moose.element( K_Pst.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = mInf/mTau*1e3
    xgate.tableB = 1/mTau*1e3

    ygate = moose.element( K_Pst.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hInf/hTau*1e3
    ygate.tableB = 1/hTau*1e3

    return K_Pst
