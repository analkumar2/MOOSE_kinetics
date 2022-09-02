# Na Persistant channel taken from mod files of Royeck2008: Nap.mod
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

def Na_P_Chan(name):
    Na = moose.HHChannel( '/library/' + name )
    Na.Ek = ENa
    Na.Gbar = 300.0*SOMA_A
    Na.Gk = 0.0
    Na.Xpower = 1.0
    Na.Ypower = 0
    Na.Zpower = 0

    # sh2   = 0
    # tha  =  -30
    # qa   = 7.2
    # Ra   = 0.4
    # Rb   = 0.124
    # thi1  = -45
    # thi2  = -45
    # qd   = 1.5
    # qg   = 1.5
    # mmin=0.02
    # hmin=0.5
    # q10=2
    # Rg   = 0.01
    # Rd   = .03
    # qq   = 10
    # tq   = -55
    # thinf  = -50
    # qinf  = 4
    # vhalfs=-60
    # a0s=0.0003
    # zetas=12
    # gms=0.2
    # smax=10
    # vvh=-58
    # vvs=2
    # a2=1
    # gbar = 0.010e4
    #
    # def trap0(v,th,a,q):
    #     if np.abs(v*1e3-th) > 1e-6:
    #         return a * (v*1e3 - th) / (1 - np.exp(-(v*1e3 - th)/q))
    #     else:
    #         return a * q
    #
    # qt=q10**((celsius-24)/10)
    # a = np.array([trap0(vm,tha+sh2,Ra,qa) for vm in v])
    # b = np.array([trap0(-vm,-tha-sh2,Rb,qa) for vm in v])
    # mtau = 1/(a+b)/qt
    # mtau[mtau<mmin] = mmin
    # minf = a/(a+b)

    minf = (1/(1+np.exp(-(v*1e3+52.3)/6.8)))
    mtau = minf*0+1


    xgate = moose.element( Na.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/mtau*1e3
    xgate.tableB = 1.0/mtau*1e3

    return Na
