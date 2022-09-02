# exec(open('../../Compilations/Kinetics/Na_Chan_Custom4.py').read())
# Na channel parameterized kinetics. Base is experimental kinetics by inf - Colbert Pan 2002. tau - Migliore2018. Plus a manual slow inacitvation gate
# Problems: q10 same for both X and Y gates.

import numpy as np
import pickle
import pandas as pd
import moose
import matplotlib.pyplot as plt

SOMA_A = 3.14e-8
F = 96485.3329
R = 8.314
celsius = 32
dt = 0.05e-3
ENa = 0.060
EK = -0.100
Eh = -0.030
ECa = 0.140
Em = -0.065

#################################
m_vhalf_inf, m_slope_inf, m_A, m_B, m_C, m_D, m_E, m_F = -31.6e-3, 6.8e-3, -3.65E-02,2.00E-02,1.61E-02,5.47E-02,3.11E-02,6.40E-04
h_vhalf_inf, h_slope_inf, h_A, h_B, h_C, h_D, h_E, h_F = -66e-3, -5.3e-3, -0.04560699,0.00433522,0.01197575,0.02617791,0.00853832,0.03900321
s_vhalf_inf, s_slope_inf, s_A, s_B, s_C, s_D, s_E, s_F = -45e-3, -6e-3, 1,0.001,0.001,0.500,0.001,1
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
    Na.Xpower = 3.0
    Na.Ypower = 1
    Na.Zpower = 1
    Na.useConcentration = 0

    [mInf,mTau] = ChanGate(v,*[m_vhalf_inf, m_slope_inf, m_A, m_B, m_C, m_D, m_E, m_F])
    [hInf,hTau] = ChanGate(v,*[h_vhalf_inf, h_slope_inf, h_A, h_B, h_C, h_D, h_E, h_F])
    [sInf,sTau] = ChanGate(v,*[s_vhalf_inf, s_slope_inf, s_A, s_B, s_C, s_D, s_E, s_F])

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

    zgate = moose.element( Na.path + '/gateZ' )
    zgate.min = Vmin
    zgate.max = Vmax
    zgate.divs = Vdivs
    zgate.tableA = sInf/sTau
    zgate.tableB = 1.0/sTau

    return Na

if __name__ == "__main__":
    [mInf,mTau] = ChanGate(v,*[m_vhalf_inf, m_slope_inf, m_A, m_B, m_C, m_D, m_E, m_F])
    [hInf,hTau] = ChanGate(v,*[h_vhalf_inf, h_slope_inf, h_A, h_B, h_C, h_D, h_E, h_F])
    # [sInf,sTau] = ChanGate(v,*[s_vhalf_inf, s_slope_inf, s_A, s_B, s_C, s_D, s_E, s_F])
    plt.figure()
    plt.plot(v, mInf, label='mInf')
    plt.plot(v, hInf, label='hInf')
    # plt.plot(v, sInf, label='sInf')
    plt.ylabel('Inf')
    plt.legend()
    plt.figure()
    plt.plot(v, mTau, label='mTau')
    plt.plot(v, hTau, label='hTau')
    # plt.plot(v, sTau, label='sTau')
    plt.ylabel('Tau')
    plt.legend()
    plt.show()

