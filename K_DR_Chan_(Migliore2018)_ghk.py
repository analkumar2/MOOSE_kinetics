# KDR channel taken from mod files of Migliore2018: kdrca1.mod
# Uses ghk instead of quasi-ohmic
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
ENa = 0.092
EK = -0.099
Eh = -0.030
ECa = 0.140
Em = -0.065

Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
# dV = (Vmax-Vmin)/Vdivs
# v = np.arange(Vmin,Vmax, dV)
v = np.linspace(Vmin,Vmax, Vdivs)
Camin = 0.04e-3
Camax = 1
Cadivs = 8000
# dCa = (Camax-Camin)/Cadivs
# ca = np.arange(Camin,Camax, dCa)
ca = np.linspace(Camin,Camax, Cadivs)

def K_DR_Chan(name):
    K_DR = moose.HHChannel( '/library/' + name )
    K_DR.Ek = EK
    K_DR.Gbar = 300.0*SOMA_A
    K_DR.Gk = 0.0
    K_DR.Xpower = 1.0
    K_DR.Ypower = 0.0
    K_DR.Zpower = 1.0
    K_DR.instant = 4

    vhalfn=13
    a0n=0.02
    zetan=-3
    gmn=0.7
    nmax=2
    q10=1
    gkdrbar=.003e4

    qt=q10**((celsius-24)/10)
    a = np.exp(1.e-3*zetan*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    ninf = 1/(1+a)
    taun = np.exp(1.e-3*zetan*gmn*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))/(qt*a0n*(1+a))
    taun[taun<nmax] = nmax/qt

    T = 33+273.15
    z = 1
    R = 8.314
    F = 96485
    Ko = 140
    DrF = z*z*F*F*Ko/R/T*v*(np.exp(z*F/R/T*(v-EK))-1)/(np.exp(z*F*v/R/T)-1)/(v-EK)

    xgate = moose.element( K_DR.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = ninf/taun*1e3
    xgate.tableB = 1.0/taun*1e3

    zgate = moose.element( K_DR.path + '/gateZ' )
    zgate.min = Vmin
    zgate.max = Vmax
    zgate.divs = Vdivs
    zgate.tableA = DrF
    zgate.tableB = DrF*0+1

    return K_DR
