# h channel taken from mod files of Migliore2018: h.mod
# Problems:

import numpy as np
import pickle
import pandas as pd
import moose

SOMA_A = 3.14e-8
F = 96485.3329
R = 8.314
celsius = 32
dt = 0.05e-3
ENa = 0.092 #from Deepanjali data
EK = -0.099 #from Deepanjali data
Eh = -0.030
ECa = 0.140 #from Deepanjali data
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

def h_Chan(name):
    h = moose.HHChannel( '/library/' + name )
    h.Ek = Eh
    h.Gbar = 300.0*SOMA_A
    h.Gk = 0.0
    h.Xpower = 1.0
    h.Ypower = 0
    h.Zpower = 0

    vhalfl=-81
    kl=-8
    vhalft=-75
    a0t=0.011
    zetat=2.2
    gmt=.4
    q10=4.5
    qtl=1
    ghdbar=.0001e4

    qt=q10**((celsius-33)/10)
    alpt = np.exp(0.0378*zetat*(v*1e3-vhalft))
    bett = np.exp(0.0378*zetat*gmt*(v*1e3-vhalft))
    a = alpt
    linf = 1/(1 + np.exp(-(v*1e3-vhalfl)/kl))
    taul = bett/(qtl*qt*a0t*(1+a))

    xgate = moose.element( h.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = linf/taul*1e3
    xgate.tableB = 1.0/taul*1e3

    return h
