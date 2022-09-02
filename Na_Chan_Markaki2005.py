# exec(open('../../Compilations/Kinetics/Na_Chan_Markaki2005.py').read())
# Na channel parameterized kinetics. Base is experimental kinetics by Makrki2005 hha2.mod
# Problems: q10 same for both X and Y gates.

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
EK = -0.100
Eh = -0.030
ECa = 0.140
Em = -0.065

#################################
m_vhalf_inf, m_slope_inf, m_A, m_B, m_C, m_D, m_E, m_F = -44e-3, 3e-3, -2.44E-06,0.5,0.5,0.5,0.5,50E-05
h_vhalf_inf, h_slope_inf, h_A, h_B, h_C, h_D, h_E, h_F = -49e-3, -3.5e-3, 1,0.13,0.33,0.18,0.05,0.004
#################################

Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
# dV = (Vmax-Vmin)/Vdivs
# v = np.arange(Vmin,Vmax, dV)
v = np.linspace(Vmin,Vmax, Vdivs)
Camin = 1e-12
Camax = 1
Cadivs = 400
# dCa = (Camax-Camin)/Cadivs
# ca = np.arange(Camin,Camax, dCa)
ca = np.linspace(Camin,Camax, Cadivs)

def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
    # alge model
    Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
    yl = (v-A)/-B
    yr = (v-A)/E
    Tau = (C + (1 + yl/(np.sqrt(1+yl**2)))/2) * (D + (1 + yr/(np.sqrt(1+yr**2)))/2) * F
    Tau[Tau<0.00002] = 0.00002
    return [Inf,Tau]

def Na_Chan(name):
    Na = moose.HHChannel( '/library/' + name )
    Na.Ek = ENa
    Na.Gbar = 300.0*SOMA_A
    Na.Gk = 0.0
    Na.Xpower = 2.0
    Na.Ypower = 1
    Na.Zpower = 0

    [mInf,mTau] = ChanGate(v,*[m_vhalf_inf, m_slope_inf, m_A, m_B, m_C, m_D, m_E, m_F])
    [hInf,hTau] = ChanGate(v,*[h_vhalf_inf, h_slope_inf, h_A, h_B, h_C, h_D, h_E, h_F])

    xgate = moose.element( Na.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = mInf/mTau
    xgate.tableB = 1.0/mTau

    ygate = moose.element( Na.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hInf/hTau
    ygate.tableB = 1.0/hTau

    return Na
