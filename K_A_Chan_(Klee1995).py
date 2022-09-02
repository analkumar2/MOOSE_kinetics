# KDR channel taken from mod files of Migliore2018: kdrca1.mod
# Custom. Kinetics changed for fitting
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

K_DR_actA = -114.1
K_DR_actB = 1.483
K_DR_actC = 0.051
K_DR_actD = -51.5
K_DR_actE = 1.2
K_DR_actF = -200
K_DR_actG = 3.84599317e-13

# # original
# K_DR_actA = -114.1
# K_DR_actB = 1.483
# K_DR_actC = 1.41182507e-01
# K_DR_actD = -7.98484941e+01
# K_DR_actE = 4.40570637e+00
# K_DR_actF = -1.14069277e+02
# K_DR_actG = 3.84599317e-13

def K_DR_Chan(name):
    K_DR = moose.HHChannel( '/library/' + name )
    K_DR.Ek = EK
    K_DR.Gbar = 300.0*SOMA_A
    K_DR.Gk = 0.0
    K_DR.Xpower = 1.0
    K_DR.Ypower = 0.0
    K_DR.Zpower = 0

    vhalfn=13
    a0n=0.02
    zetan=-3
    gmn=0.7
    nmax=2
    q10=1
    gkdrbar=.003e4

    # qt=q10**((celsius-24)/10)
    # a = np.exp(1.e-3*zetan*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    # ninf = 1/(1+a)
    # taun = np.exp(1.e-3*zetan*gmn*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))/(qt*a0n*(1+a))
    # # taun[taun<nmax] = nmax/qt
    # taun = taun*1e-3

    ninf = 1/(1+np.exp(K_DR_actA*v+K_DR_actB))
    # taun = 1e3/(K_DR_actC*v*np.exp(K_DR_actD*v)+K_DR_actE*np.exp(K_DR_actF*v))
    taun=K_DR_actC*np.exp(K_DR_actD*v)/(1+K_DR_actE*np.exp(K_DR_actF*v)) + K_DR_actG

    xgate = moose.element( K_DR.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = ninf/taun
    xgate.tableB = 1.0/taun

    return K_DR
