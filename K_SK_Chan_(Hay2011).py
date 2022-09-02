# K_SK channelSK_E2.mod Hay2011
# exec(open('K_SK_Chan_(Hay2011).py').read())

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

def K_SK_Chan(name):
    K_SK = moose.HHChannel( '/library/' + name )
    K_SK.Ek = EK
    K_SK.Gbar = 300.0*SOMA_A
    K_SK.Gk = 0.0
    K_SK.Xpower = 0.0
    K_SK.Ypower = 0.0
    K_SK.Zpower = 1.0
    K_SK.useConcentration = 1

    zInf = 1/(1 + (0.00043 / ca)**4.8)
    zTau = 1 + 0*ca

    zgate = moose.element( K_SK.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zgate.tableA = zInf/zTau*1e3
    zgate.tableB = 1/zTau*1e3

    addmsg3 = moose.Mstring( K_SK.path + '/addmsg3' )
    addmsg3.value = '../Ca_conc concOut . concen'

    return K_SK
