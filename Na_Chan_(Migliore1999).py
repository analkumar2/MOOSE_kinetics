# exec(open('../../Compilations/Kinetics/Na_Chan_(Migliore1999).py').read())
# Na channel taken from mod files of Migliore1999: na3n.mod
# Fitted version

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

def Na_Chan(name):
    Na = moose.HHChannel( '/library/' + name )
    Na.Ek = ENa
    Na.Gbar = 300.0*SOMA_A
    Na.Gk = 0.0
    Na.Xpower = 3.0
    Na.Ypower = 1
    Na.Zpower = 0

    def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
        # alge model
        Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
        yl = (v-A)/-B
        yr = (v-A)/E
        Tau = (C + (1 + yl/(np.sqrt(1+yl**2)))/2) * (D + (1 + yr/(np.sqrt(1+yr**2)))/2) * F
        Tau[Tau<0.00002] = 0.00002
        return [Inf,Tau]

    minf, mtau = ChanGate(v,-0.03843252,0.0072, -0.03737914,0.02071916,0.01278186,0.05709504,0.02957831,0.0038622)
    hinf, htau = ChanGate(v,-0.05,-0.004, -0.04590684,0.004695,0.00596418,0.0263036,0.0081486,0.23889591)

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
