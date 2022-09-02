# As seen in Migliore2018
# Just to setup. The actual values are set up after model building.
# Calcium pool decay and conversion from calcium current to concentration.
# set moose.element('/model/elec/soma/Ca_conc').B = 1000e3/(2*F*depth*np.pi*sm_diam*sm_len*2)

import numpy as np
import pickle
import pandas as pd
import moose

#################################
taur    = 100e-3
cainf    = 0.05e-3
#################################

def Ca_Conc ( name ):
    Ca = moose.CaConc( '/library/' + name )
    depth    = 0.1e-6
    # B = 28789637.7

    Ca.tau = taur
    Ca.Ca_base = cainf
    Ca.diameter = 20e-6
    Ca.length = 20e-6
    Ca.thick = 0.1e-6
    return Ca
