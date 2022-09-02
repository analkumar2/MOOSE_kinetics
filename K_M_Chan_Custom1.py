# K_M channel Im.mod Hay2011
# exec(open('K_M_Chan_(Hay2011).py').read())

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

############################################
factor = 3.3e-3
############################################


def K_M_Chan(name):
    K_M = moose.HHChannel( '/library/' + name )
    K_M.Ek = EK
    K_M.Gbar = 300.0*SOMA_A
    K_M.Gk = 0.0
    K_M.Xpower = 1.0
    K_M.Ypower = 0.0
    K_M.Zpower = 0.0

    qt = 2.3**((34-21)/10)
    mAlpha = factor*np.exp(2.5*0.04*(v*1e3 - -35))
    mBeta = factor*np.exp(-2.5*0.04*(v*1e3 - -35))

    xgate = moose.element( K_M.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = mAlpha*1e3
    xgate.tableB = qt*(mAlpha+mBeta)*1e3

    return K_M
