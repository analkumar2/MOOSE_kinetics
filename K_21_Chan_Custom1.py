# exec(open('../../Compilations/Kinetics/K_21_Chan_Custom1.py').read())
# Channelpedia kinetics after fitting. fit_Kv2.1_mouse_CHO35_32675

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

[minf_m90, minf_m80, minf_m70, minf_m60, minf_m50, minf_m40, minf_m30, minf_m20, minf_m10, minf_0, minf_10, minf_20, minf_30, minf_40, minf_50, minf_60, minf_70, minf_80] = [0, 0, 0, 0, 0, 0, 0.014778557902240114, 0.12440907737836038, 0.32644969164762005, 0.4939673946612481, 0.6062082677794262, 0.6925405431544853, 0.7547172710411801, 0.8041219928637904, 0.8357290533410773, 0.8624147007018226, 0.8841496013479486, 0.8979634200234778]
[mtau_m90, mtau_m80, mtau_m70, mtau_m60, mtau_m50, mtau_m40, mtau_m30, mtau_m20, mtau_m10, mtau_0, mtau_10, mtau_20, mtau_30, mtau_40, mtau_50, mtau_60, mtau_70, mtau_80] = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.09999999996162072, 0.020777818330050365, 0.010671060904288474, 0.006703769145563319, 0.0044309113387414, 0.003224061657330687, 0.002524784147194683, 0.00210592809771354, 0.0018592671524697525, 0.0016749012882572218, 0.0015627450318380699, 0.0014749072662173194]
[hinf_m90, hinf_m80, hinf_m70, hinf_m60, hinf_m50, hinf_m40, hinf_m30, hinf_m20, hinf_m10, hinf_0, hinf_10, hinf_20, hinf_30, hinf_40, hinf_50, hinf_60, hinf_70, hinf_80] = [1, 1, 1, 1, 1, 0.3, 2.253655822247467e-10, 2.3695543557633294e-17, 0.05264001822205584, 0.13336298143641784, 0.1474047523692984, 0.1647609967027383, 0.19287397670410475, 0.21015558524246095, 0.22657341896934804, 0.23777997846836482, 0.24949074522445402, 0.2563905322354998]
[htau_m90, htau_m80, htau_m70, htau_m60, htau_m50, htau_m40, htau_m30, htau_m20, htau_m10, htau_0, htau_10, htau_20, htau_30, htau_40, htau_50, htau_60, htau_70, htau_80] = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2775366296151369, 0.3972217370547666, 0.2748642406543788, 0.20193562149025931, 0.18691986665827, 0.17758044785985413, 0.16568374298787542, 0.16122751315663916, 0.1572864085873158, 0.15557833831912238, 0.1528, 0.15]
Xpower, Ypower = 1,1

def K_21_Chan(name):
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