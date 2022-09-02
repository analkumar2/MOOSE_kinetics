# exec(open('../../Compilations/Kinetics/K_42_Chan_Custom1.py').read())
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

[minf_m90, minf_m80, minf_m70, minf_m60, minf_m50, minf_m40, minf_m30, minf_m20, minf_m10, minf_0, minf_10, minf_20, minf_30, minf_40, minf_50, minf_60, minf_70, minf_80] = [0, 0, 0, 0, 0, 0, 5.97265239177511e-14, 3.328565161997558e-11, 0.3, 0.5429664374500748, 0.5627165859088421, 0.5559579996396777, 0.566332215396569, 0.5142681553571444, 0.49732164093782705, 0.012185281688337815, 0.1837140070988242, 0.024128789742599216]
[mtau_m90, mtau_m80, mtau_m70, mtau_m60, mtau_m50, mtau_m40, mtau_m30, mtau_m20, mtau_m10, mtau_0, mtau_10, mtau_20, mtau_30, mtau_40, mtau_50, mtau_60, mtau_70, mtau_80] = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 4.1308665584688375e-06, 0.00015219051506553568, 0.01, 0.00039010418946787987, 0.00034341588838282766, 0.00024030835772276895, 0.0001774694081040018, 0.0001198632926717059, 8.734656017332724e-05, 0.09073677454608732, 2.347339858711643e-09, 9.040274651119406e-06]
[hinf_m90, hinf_m80, hinf_m70, hinf_m60, hinf_m50, hinf_m40, hinf_m30, hinf_m20, hinf_m10, hinf_0, hinf_10, hinf_20, hinf_30, hinf_40, hinf_50, hinf_60, hinf_70, hinf_80] = [1, 1, 1, 1, 1, 0.3, 7.189764409679599e-12, 3.227968267336322e-05, 0.3, 6.168750596444721e-14, 8.779505665440972e-13, 0.005545246367651815, 0.0042715646152740585, 0.009740701435088902, 0.016980739899913454, 0.9600066036748729, 0.07757601951126956, 0.09207713436782324]
[htau_m90, htau_m80, htau_m70, htau_m60, htau_m50, htau_m40, htau_m30, htau_m20, htau_m10, htau_0, htau_10, htau_20, htau_30, htau_40, htau_50, htau_60, htau_70, htau_80] = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.820491652392009e-11, 0.036741265987020405, 0.1, 0.009247638013742109, 0.006426495363261932, 0.005513542614026591, 0.0049219833223978, 0.004802004051211955, 0.0046118665724215096, 6.494521914075398, 0.004410634822090242, 2.8306227788350165]
Xpower, Ypower = 1,1

def K_42_Chan(name):
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