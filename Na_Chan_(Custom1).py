# Na channel custom
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
Camin = 1e-12
Camax = 1
Cadivs = 400
# dCa = (Camax-Camin)/Cadivs
# ca = np.arange(Camin,Camax, dCa)
ca = np.linspace(Camin,Camax, Cadivs)

Na_actA = -138.9
Na_actB = -5.34
Na_actC = 2.60187360e-05
Na_actD = -7.05694169e+01
Na_actE = 2.64467199e-02
Na_actF = -1.02623323e+02
Na_actG = 1.73923802e-05

Na_inactA = 350
Na_inactB = 12.5
Na_inactC = 2.004014672869906e-09
Na_inactD = -360.599
Na_inactE = 1.866086164497956e-09
Na_inactF = -454.5
Na_inactG = 0.00047

# # original
# Na_actA = -138.9
# Na_actB = -5.34
# Na_actC = 2.60187360e-05
# Na_actD = -7.05694169e+01
# Na_actE = 2.64467199e-02
# Na_actF = -1.02623323e+02
# Na_actG = 1.73923802e-05
#
# Na_inactA = 250
# Na_inactB = 12.5
# Na_inactC = 2.004014672869906e-09
# Na_inactD = -360.599
# Na_inactE = 1.866086164497956e-09
# Na_inactF = -454.5
# Na_inactG = 0.00047

def Na_Chan(name):
    Na = moose.HHChannel( '/library/' + name )
    Na.Ek = ENa
    Na.Gbar = 300.0*SOMA_A
    Na.Gk = 0.0
    Na.Xpower = 3.0
    Na.Ypower = 1
    Na.Zpower = 0

    sh2   = 0
    tha  =  -30
    qa   = 7.2
    Ra   = 0.4
    Rb   = 0.124
    thi1  = -45
    thi2  = -45
    qd   = 1.5
    qg   = 1.5
    mmin=0.02
    hmin=0.5
    q10=2
    Rg   = 0.01
    Rd   = .03
    qq   = 10
    tq   = -55
    thinf  = -50
    qinf  = 4
    vhalfs=-60
    a0s=0.0003
    zetas=12
    gms=0.2
    smax=10
    vvh=-58
    vvs=2
    a2=1
    gbar = 0.010e4

    def trap0(v,th,a,q):
        if np.abs(v*1e3-th) > 1e-6:
            return a * (v*1e3 - th) / (1 - np.exp(-(v*1e3 - th)/q))
        else:
            return a * q

    # qt=q10**((celsius-24)/10)
    # a = np.array([trap0(vm,tha+sh2,Ra,qa) for vm in v])
    # b = np.array([trap0(-vm,-tha-sh2,Rb,qa) for vm in v])
    # mtau = 1/(a+b)/qt
    # # mtau[mtau<mmin] = mmin
    # minf = a/(a+b)
    # mtau = mtau*1e-3
    #
    # a = np.array([trap0(vm,thi1+sh2,Rd,qd) for vm in v])
    # b = np.array([trap0(-vm,-thi2-sh2,Rg,qg) for vm in v])
    # htau =  1/(a+b)/qt
    # # htau[htau<hmin] = hmin
    # hinf = 1/(1+np.exp((v*1e3-thinf-sh2)/qinf))
    # htau = htau*1e-3

    minf = 1/(1+np.exp(Na_actA*v+Na_actB))
    # taun = 1e3/(Na_actC*v*np.exp(Na_actD*v)+Na_actE*np.exp(Na_actF*v))
    mtau=Na_actC*np.exp(Na_actD*v)/(1+Na_actE*np.exp(Na_actF*v)) + Na_actG

    hinf = 1/(1+np.exp(Na_inactA*v+Na_inactB))
    # taun = 1e3/(Na_inactC*v*np.exp(Na_inactD*v)+Na_inactE*np.exp(Na_inactF*v))
    htau=Na_inactC*np.exp(Na_inactD*v)/(1+Na_inactE*np.exp(Na_inactF*v)) + Na_inactG


    xgate = moose.element( Na.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/mtau
    xgate.tableB = 1.0/mtau

    ygate = moose.element( Na.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/htau
    ygate.tableB = 1.0/htau

    return Na
