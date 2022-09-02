# KSK channel taken from mod files of Migliore2018: kca.mod
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
Camin = 1e-12 #If changing this, be careful that the dCa is at most 0.01e-3
Camax = 3e-3 #3e-3 works because the tables becomes steady after 3e-3 for K_SK
Cadivs = 4000
# dCa = (Camax-Camin)/Cadivs
# ca = np.arange(Camin,Camax, dCa)
ca = np.linspace(Camin,Camax, Cadivs)

def K_SK_Chan(name):
    K_SK = moose.HHChannel( '/library/' + name )
    K_SK.Ek = EK
    K_SK.Gbar = 300.0*SOMA_A
    K_SK.Gk = 0.0
    K_SK.Xpower = 0
    K_SK.Ypower = 0
    K_SK.Zpower = 3
    K_SK.useConcentration = 1

    beta    = 0.03
    cac     = 0.00035
    taumin  = 0.5
    q10 = 3
    gbar    = 0.01e4

    tadj = q10**((celsius-22.0)/10)
    car = (ca/cac)**4
    m_inf = car / ( 1 + car )
    tau_m =  1 / beta / (1 + car) / tadj
    tau_m[tau_m<taumin] = taumin

    zgate = moose.element( K_SK.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zgate.tableA = m_inf/tau_m*1e3
    zgate.tableB = 1.0/tau_m*1e3

    addmsg3 = moose.Mstring( K_SK.path + '/addmsg3' )
    addmsg3.value = '../Ca_conc concOut . concen'
    return K_SK
