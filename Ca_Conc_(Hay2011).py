# As seen in Hay2011
# Just to setup. The actual values are set up after model building.
# Calcium pool decay and conversion from calcium current to concentration.
# set moose.element('/model/elec/soma/Ca_conc').B = 10000*gamma/(2*F*depth)*1e4/(sm_area) where default gamma=0.05 and depth=0.1

import numpy as np
import pickle
import pandas as pd
import moose

def Ca_Conc ( name ):
    Ca = moose.CaConc( '/library/' + name )
    depth    = 0.1e-6 #Doesnt' matter
    taur    = 80e-3
    cainf    = 0.1e-3 #Higher than the standard 0.05e-3
    #B = 252/sm_area

    Ca.tau = taur
    Ca.Ca_base = cainf
    Ca.diameter = 20e-6 #Doesn't matter
    Ca.length = 20e-6 #Doesn't matter
    Ca.thick = 0.1e-6 #Doesn't matter
    return Ca
