# HVA type calcium channel Ca_HVA.mod Hay2011
# exec(open('Ca_HVA_Chan_(Hay2011).py').read())

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

def Ca_HVA_Chan(name):
    Ca_HVA = moose.HHChannel( '/library/' + name )
    Ca_HVA.Ek = ECa
    Ca_HVA.Gbar = 300.0*SOMA_A
    Ca_HVA.Gk = 0.0
    Ca_HVA.Xpower = 2.0
    Ca_HVA.Ypower = 1.0
    Ca_HVA.Zpower = 0.0

    mAlpha =  (0.055*(-27-v*1e3))/(np.exp((-27-v*1e3)/3.8) - 1)
    mBeta  =  (0.94*np.exp((-75-v*1e3)/17))
    hAlpha =  (0.000457*np.exp((-13-v*1e3)/50))
    hBeta  =  (0.0065/(np.exp((-v*1e3-15)/28)+1))

    xgate = moose.element( Ca_HVA.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = mAlpha*1e3
    xgate.tableB = (mAlpha + mBeta)*1e3

    ygate = moose.element( Ca_HVA.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hAlpha*1e3
    ygate.tableB = (hAlpha + hBeta)*1e3

    addmsg2 = moose.Mstring( Ca_HVA.path + '/addmsg2' )
    addmsg2.value = '. IkOut ../Ca_conc current'
    return Ca_HVA
