# L type calcium channel cal2.mod Migliore2018
# exec(open('Ca_L_Chan_Custom1.py').read())

import numpy as np
import pickle
import pandas as pd
import moose
import sys

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

###############################
a0m = 0.1
###############################

Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
# dV = (Vmax-Vmin)/Vdivs
# v = np.arange(Vmin,Vmax, dV)
v = np.linspace(Vmin,Vmax, Vdivs)
Camin = 1e-12
Camax = 1
Cadivs = 10000 # Is enough for Ca_L
# dCa = (Camax-Camin)/Cadivs
# ca = np.arange(Camin,Camax, dCa)
ca = np.linspace(Camin,Camax, Cadivs)

def Ca_L_Chan(name):
    Ca_L = moose.HHChannel2D( '/library/' + name )
    Ca_L.Ek = ECa
    Ca_L.Gbar = 300.0*SOMA_A
    Ca_L.Gk = 0.0
    Ca_L.Xpower = 2.0
    Ca_L.Ypower = 0.0
    Ca_L.Zpower = 1.0
    Ca_L.Xindex = 'VOLT_INDEX'
    # Ca_L.Yindex = 'C1_INDEX'
    Ca_L.Zindex = 'VOLT_C1_INDEX'
    Ca_L.instant = 4

    xgate = moose.element( Ca_L.path + '/gateX' )
    xgate.xminA = Vmin
    xgate.xmaxA = Vmax
    xgate.xdivsA = Vdivs
    # xgate.yminA = Camin
    # xgate.ymaxA = Camax
    # xgate.ydivsA = Cadivs
    xgate.xminB = Vmin
    xgate.xmaxB = Vmax
    xgate.xdivsB = Vdivs
    # xgate.yminB = Camin
    # xgate.ymaxB = Camax
    # xgate.ydivsB = Cadivs

    zgate = moose.element( Ca_L.path + '/gateZ' )
    zgate.xminA = Vmin
    zgate.xmaxA = Vmax
    zgate.xdivsA = Vdivs
    zgate.yminA = Camin
    zgate.ymaxA = Camax
    zgate.ydivsA = Cadivs
    zgate.xminB = Vmin
    zgate.xmaxB = Vmax
    zgate.xdivsB = Vdivs
    zgate.yminB = Camin
    zgate.ymaxB = Camax
    zgate.ydivsB = Cadivs

    cao = 2
    z = 2
    T = celsius+273.15
    ki=.001
    cai = 50.e-6
    cao = 2
    q10 = 5
    mmin=0.2
    tfa = 1
    # a0m =0.1
    zetam = 2
    vhalfm = 4
    gmm=0.1
    gcalbar=.003e4

    qt=q10**((celsius-25)/10)
    a = 15.69*(-1.0*v*1e3+81.5)/(np.exp((-1.0*v*1e3+81.5)/10.0)-1.0)
    b = 1/(a + 0.29*np.exp(-v*1e3/10.86))
    minf = a*b
    alpmt = np.exp(0.0378*zetam*(v*1e3-vhalfm))
    betmt = np.exp(0.0378*zetam*gmm*(v*1e3-vhalfm))
    tau = betmt/(qt*a0m*(1+alpmt))
    tau[tau<mmin/qt]=mmin/qt
    XtblA = np.zeros([Vdivs,1])
    XtblB = np.zeros([Vdivs,1])
    for i in np.arange(1):
        XtblA[:,i] = minf/tau
        XtblB[:,i] = 1/tau
        print(i, end='\r')
    xgate.tableA = XtblA*1e3
    xgate.tableB = XtblB*1e3

    # ezfrt  = np.exp(z*v*F/R/T)
    # ZtblA = np.zeros([Vdivs,Cadivs])
    # ZtblB = np.zeros([Vdivs,Cadivs])
    # for i in np.arange(Cadivs):
    #     ZtblA[:,i] = ki/(ki+ca[i])*v*(ca[i]/cao*ezfrt-1)/(ezfrt-1)/(v-ECa)
    #     print(i, end='\r')
    # ZtblB = ZtblA*0 +1
    # zgate.tableA = ZtblA
    # zgate.tableB = ZtblB

    ezfrt  = np.exp(z*v*F/R/T)
    caezfrt = np.outer(ca,ezfrt)/cao - 1
    caezfrt_ex = v/(ezfrt-1)/(v-ECa)
    ZtblA = np.transpose(caezfrt*caezfrt_ex)*(ki/(ki+ca))
    ZtblB = ZtblA*0 +1
    zgate.tableA = ZtblA
    zgate.tableB = ZtblB

    # np.savetxt('temp.txt',np.isclose(ZtblA, ZtblA_new))

    # np.savez('Ca_L_Chan_Custom1_tbls.npz', XtblA=XtblA, XtblB=XtblB, ZtblA=ZtblA.ravel(), ZtblB=ZtblB.ravel())
    # tbls = np.load('Ca_L_Chan_Custom1_tbls.npz')
    # xgate.tableA = tbls['XtblA']
    # xgate.tableB = tbls['XtblB']
    # zgate.tableA = tbls['ZtblA'].reshape(Vdivs,Cadivs)
    # zgate.tableB = tbls['ZtblB'].reshape(Vdivs,Cadivs)

    addmsg4 = moose.Mstring( Ca_L.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc concOut . concen'

    addmsg2 = moose.Mstring( Ca_L.path + '/addmsg2' )
    addmsg2.value = '. IkOut ../Ca_conc current'
    return Ca_L
