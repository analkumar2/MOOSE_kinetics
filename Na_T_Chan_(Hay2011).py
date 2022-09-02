# Na_T channel NaTs2_t.mod Hay2011
# exec(open('Na_T_Chan_(Hay2011).py').read())

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

def Na_T_Chan(name):
    Na_T = moose.HHChannel( '/library/' + name )
    Na_T.Ek = ENa
    Na_T.Gbar = 300.0*SOMA_A
    Na_T.Gk = 0.0
    Na_T.Xpower = 3.0
    Na_T.Ypower = 1.0
    Na_T.Zpower = 0.0

    qt = 2.3**((34-21)/10)
    mAlpha = (0.182 * (v*1e3- -32))/(1-(np.exp(-(v*1e3- -32)/6)))
    mBeta  = (0.124 * (-v*1e3 -32))/(1-(np.exp(-(-v*1e3 -32)/6)))
    mInf = mAlpha/(mAlpha + mBeta)
    mTau = (1/(mAlpha + mBeta))/qt
    hAlpha = (-0.015 * (v*1e3- -60))/(1-(np.exp((v*1e3- -60)/6)))
    hBeta  = (-0.015 * (-v*1e3 -60))/(1-(np.exp((-v*1e3 -60)/6)))
    hInf = hAlpha/(hAlpha + hBeta)
    hTau = (1/(hAlpha + hBeta))/qt

    xgate = moose.element( Na_T.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = mInf/mTau*1e3
    xgate.tableB = 1/mTau*1e3

    ygate = moose.element( Na_T.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hInf/hTau*1e3
    ygate.tableB = 1/hTau*1e3

    return Na_T
