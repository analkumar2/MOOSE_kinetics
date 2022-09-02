# exec(open('../../Compilations/Kinetics/Na_Chan_(Poolos2002).py').read())
# Na channel taken from mod files of Poolos2002: na3n.mod
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

    minf, mtau = ChanGate(v,-0.02343252,0.0072, -0.02148724,0.02037766,0.01387742,0.05055026,0.03167365,0.00384101)
    hinf, htau = ChanGate(v,-0.035,-0.004, -0.03067039,0.00462528,0.00659626,0.02280288,0.0087322,0.23511362)

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
