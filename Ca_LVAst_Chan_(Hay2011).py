# HVA type calcium channel Ca_LVAst.mod Hay2011
# exec(open('Ca_LVAst_Chan_(Hay2011).py').read())

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

def Ca_LVAst_Chan(name):
    Ca_LVAst = moose.HHChannel( '/library/' + name )
    Ca_LVAst.Ek = ECa
    Ca_LVAst.Gbar = 300.0*SOMA_A
    Ca_LVAst.Gk = 0.0
    Ca_LVAst.Xpower = 2.0
    Ca_LVAst.Ypower = 1.0
    Ca_LVAst.Zpower = 0.0

    qt = 2.3**((34-21)/10)
    V = v+0.010 #Because Hay2011 shifted this
    mInf = 1.0000/(1+ np.exp((V*1e3 - -30.000)/-6))
    mTau = (5.0000 + 20.0000/(1+np.exp((V*1e3 - -25.000)/5)))/qt
    hInf = 1.0000/(1+ np.exp((V*1e3 - -80.000)/6.4))
    hTau = (20.0000 + 50.0000/(1+np.exp((V*1e3 - -40.000)/7)))/qt

    xgate = moose.element( Ca_LVAst.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = mInf/mTau*1e3
    xgate.tableB = 1/mTau*1e3

    ygate = moose.element( Ca_LVAst.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hInf/hTau*1e3
    ygate.tableB = 1/hTau*1e3

    addmsg2 = moose.Mstring( Ca_LVAst.path + '/addmsg2' )
    addmsg2.value = '. IkOut ../Ca_conc current'
    return Ca_LVAst
