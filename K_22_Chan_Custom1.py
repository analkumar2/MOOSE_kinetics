# exec(open('../../Compilations/Kinetics/K_22_Chan_Custom1.py').read())
# Channelpedia kinetics after fitting. fit_Kv2.2_rat_CHO35_7666

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

SOMA_A = 3.14e-8
F = 96485.3329
R = 8.314
celsius = 32
dt = 0.05e-3

Vmin = -0.090
Vmax = 0.080
Vdivs = 18
v = np.linspace(Vmin,Vmax, Vdivs)

[minf_m90, minf_m80, minf_m70, minf_m60, minf_m50, minf_m40, minf_m30, minf_m20, minf_m10, minf_0, minf_10, minf_20, minf_30, minf_40, minf_50, minf_60, minf_70, minf_80] = [0, 0, 0, 0, 0, 0, 0.14537763030775508, 0.3098279802117982, 0.4480877075658604, 0.5569458932103868, 0.6268796006736828, 0.6881757867589934, 0.7391385069304409, 0.7833163065281092, 0.8294264707631671, 0.8604407457240514, 0.8892205686252913, 0.9105871408756966]
[mtau_m90, mtau_m80, mtau_m70, mtau_m60, mtau_m50, mtau_m40, mtau_m30, mtau_m20, mtau_m10, mtau_0, mtau_10, mtau_20, mtau_30, mtau_40, mtau_50, mtau_60, mtau_70, mtau_80] = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.01045701930219431, 0.006442508082027577, 0.003960080931819326, 0.002936002161377107, 0.002300537523748545, 0.0018509333454103667, 0.0015540985890712406, 0.0013634465175197225, 0.0012238990098193175, 0.0011328242697070216, 0.0010569961470974737, 0.0009632584441008515]
[hinf_m90, hinf_m80, hinf_m70, hinf_m60, hinf_m50, hinf_m40, hinf_m30, hinf_m20, hinf_m10, hinf_0, hinf_10, hinf_20, hinf_30, hinf_40, hinf_50, hinf_60, hinf_70, hinf_80] = [1, 1, 1, 1, 1, 0.3, 0.03790161435438773, 0.08990959911917082, 0.13495284064299926, 0.18577329497978826, 0.2132079836951254, 0.23344455686133914, 0.2432487364859475, 0.25344805294216066, 0.2689671129407693, 0.2895260663642892, 0.30719387467604753, 0.31774394900458086]
[htau_m90, htau_m80, htau_m70, htau_m60, htau_m50, htau_m40, htau_m30, htau_m20, htau_m10, htau_0, htau_10, htau_20, htau_30, htau_40, htau_50, htau_60, htau_70, htau_80] = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.41927884508025026, 0.27838983174667487, 0.23349884137457239, 0.20475178442318262, 0.2066139009982512, 0.20969008422316998, 0.2144315953207104, 0.21850766994905096, 0.20908802055689718, 0.20309744036745572, 0.19595635032873898, 0.19560528797087373]
Xpower, Ypower = 1,1

def K_22_Chan(name):
    if moose.exists('/library/' + name):
        K = moose.element('/library/' + name)
    else:
        K = moose.HHChannel( '/library/' + name )
    K.Ek = EK
    K.Gbar = 300.0*SOMA_A
    K.Gk = 0.0
    K.Xpower = Xpower
    K.Ypower = Ypower
    K.Zpower = 0
    K.useConcentration = 0

    mInf = np.array([minf_m90, minf_m80, minf_m70, minf_m60, minf_m50, minf_m40, minf_m30, minf_m20, minf_m10, minf_0, minf_10, minf_20, minf_30, minf_40, minf_50, minf_60, minf_70, minf_80])
    mTau = np.array([mtau_m90, mtau_m80, mtau_m70, mtau_m60, mtau_m50, mtau_m40, mtau_m30, mtau_m20, mtau_m10, mtau_0, mtau_10, mtau_20, mtau_30, mtau_40, mtau_50, mtau_60, mtau_70, mtau_80])
    hInf = np.array([hinf_m90, hinf_m80, hinf_m70, hinf_m60, hinf_m50, hinf_m40, hinf_m30, hinf_m20, hinf_m10, hinf_0, hinf_10, hinf_20, hinf_30, hinf_40, hinf_50, hinf_60, hinf_70, hinf_80])
    hTau = np.array([htau_m90, htau_m80, htau_m70, htau_m60, htau_m50, htau_m40, htau_m30, htau_m20, htau_m10, htau_0, htau_10, htau_20, htau_30, htau_40, htau_50, htau_60, htau_70, htau_80])

    xgate = moose.element( K.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = mInf/mTau
    xgate.tableB = 1.0/mTau

    ygate = moose.element( K.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hInf/hTau
    ygate.tableB = 1.0/hTau

    return K


if __name__ == "__main__":
    nInf = [minf_m90, minf_m80, minf_m70, minf_m60, minf_m50, minf_m40, minf_m30, minf_m20, minf_m10, minf_0, minf_10, minf_20, minf_30, minf_40, minf_50, minf_60, minf_70, minf_80]
    nTau = [mtau_m90, mtau_m80, mtau_m70, mtau_m60, mtau_m50, mtau_m40, mtau_m30, mtau_m20, mtau_m10, mtau_0, mtau_10, mtau_20, mtau_30, mtau_40, mtau_50, mtau_60, mtau_70, mtau_80]
    lInf = [hinf_m90, hinf_m80, hinf_m70, hinf_m60, hinf_m50, hinf_m40, hinf_m30, hinf_m20, hinf_m10, hinf_0, hinf_10, hinf_20, hinf_30, hinf_40, hinf_50, hinf_60, hinf_70, hinf_80]
    lTau = [htau_m90, htau_m80, htau_m70, htau_m60, htau_m50, htau_m40, htau_m30, htau_m20, htau_m10, htau_0, htau_10, htau_20, htau_30, htau_40, htau_50, htau_60, htau_70, htau_80]

    plt.figure()
    plt.plot(v, nInf, label='nInf')
    plt.plot(v, lInf, label='lInf')
    plt.ylabel('Inf')
    plt.legend()
    plt.figure()
    plt.plot(v, nTau, label='nTau')
    plt.plot(v, lTau, label='lTau')
    plt.ylabel('Tau')
    plt.legend()
    plt.show()