# exec(open('../../Compilations/Kinetics/K_A_Chan_Custom3.py').read())
# Migliore2018 kinetics. Parameterized but inactivation tau not by alge model

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
ENa = 0.092 #from Deepanjali data
EK = -0.099 #from Deepanjali data
Eh = -0.030
ECa = 0.140 #from Deepanjali data
Em = -0.065


#################################
#Not corrected for junction potential (-7e-3V). Temp was 35*C.
n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F = 0.0112,0.017, -8.78e-3,5.63e-2,0,0,2.65e-2,1.05e-2
l_vhalf_inf, l_slope_inf, l_min, l_m, l_cm = -0.056,-0.00877, 0.002, 0.26, 0.050
#################################

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

def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
    # alge model
    Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
    yl = (v-A)/-B
    yr = (v-A)/E
    Tau = (C + (1 + yl/(np.sqrt(1+yl**2)))/2) * (D + (1 + yr/(np.sqrt(1+yr**2)))/2) * F
    Tau[Tau<0.00002] = 0.00002
    return [Inf,Tau]

def K_A_Chan(name):
    K_A = moose.HHChannel( '/library/' + name )
    K_A.Ek = EK
    K_A.Gbar = 300.0*SOMA_A
    K_A.Gk = 0.0
    K_A.Xpower = 1.0
    K_A.Ypower = 1.0
    K_A.Zpower = 0

    [nInf,nTau] = ChanGate(v,*[n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])
    lInf = ChanGate(v,*[l_vhalf_inf, l_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])[0]


    lTau = l_m*(v+l_cm)
    lTau[lTau<l_min] = l_min


    xgate = moose.element( K_A.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = nInf/nTau
    xgate.tableB = 1.0/nTau

    ygate = moose.element( K_A.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = lInf/lTau
    ygate.tableB = 1.0/lTau

    return K_A


if __name__ == "__main__":
    [nInf,nTau] = ChanGate(v,*[n_vhalf_inf, n_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])
    lInf = ChanGate(v,*[l_vhalf_inf, l_slope_inf, n_A, n_B, n_C, n_D, n_E, n_F])[0]
    lTau = l_m*(v+l_cm)
    lTau[lTau<l_min] = l_min

    # hhh = np.arange(-55e-3 - 0.016, 36e-3 - 0.016, 10e-3)
    # parameterfitted_list = np.array([np.array([2.56015621e-09, 1.12713826e-02, 3.82738005e-18, 5.00000000e-01,
    #    2.43722849e-09, 2.81531338e-01]), np.array([5.01958073e-09, 1.30315203e-02, 8.51866701e-04, 1.14235667e-01,
    #    4.33268931e-09, 1.14243302e-01]), np.array([7.37120030e-09, 1.25073616e-02, 1.36142786e-06, 9.74136065e-02,
    #    5.34766250e-09, 9.74163312e-02]), np.array([8.99721426e-09, 9.24950966e-03, 5.35874431e-03, 1.49499585e-01,
    #    6.20488857e-09, 1.49505473e-01]), np.array([1.16121395e-08, 7.32873273e-03, 7.26493010e-11, 3.46096004e-01,
    #    8.67817670e-09, 5.00000000e-01]), np.array([1.54864626e-08, 7.78551206e-03, 7.15086937e-10, 8.17741596e-02,
    #    1.29630899e-08, 7.50310292e-02]), np.array([1.24100284e-08, 1.57381849e-01, 6.50940375e-02, 5.36024556e-02,
    #    1.76802248e-08, 6.50738526e-03]), np.array([1.53142103e-08, 5.27067710e-03, 9.55547838e-01, 4.99999411e-01,
    #    7.90179057e-09, 1.17767958e-02]), np.array([2.54855042e-08, 6.10498775e-03, 9.12926651e-01, 5.00000000e-01,
    #    2.17731422e-09, 6.85336947e-02]), np.array([2.92396179e-08, 6.27161627e-03, 9.02253739e-01, 5.00000000e-01,
    #    1.56575057e-09, 1.07168695e-01])])
    # nInf_fitted = parameterfitted_list[:,0]/3e-8
    # nTau_fitted = parameterfitted_list[:,1]

    # lInf_fitted = parameterfitted_list[:,2]
    # lTau_fitted = parameterfitted_list[:,3]

    plt.figure()
    plt.plot(v, nInf, label='nInf')
    plt.plot(v, lInf, label='lInf')
    # plt.plot(hhh, nInf_fitted, label='nInf_fitted')
    # plt.plot(hhh, lInf_fitted, label='lInf_fitted')
    plt.ylabel('Inf')
    plt.legend()
    plt.figure()
    plt.plot(v, nTau, label='nTau')
    plt.plot(v, lTau, label='lTau')
    # plt.plot(hhh, nTau_fitted, label='nTau_fitted')
    # plt.plot(hhh, lTau_fitted, label='lTau_fitted')
    plt.ylabel('Tau')
    plt.legend()
    plt.show()