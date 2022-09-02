import moose
import numpy as np
import matplotlib.pyplot as plt
from Ca_L_Chan_Custom1 import Ca_L_Chan
moose.Neutral('/library')
Ca_L_Chan('Ca_L_Chan')
#tbls = np.load('Ca_L_Chan_Custom1_tbls.npz')
# tableA = tbls['XtblA']
# tableB = tbls['XtblB']
# #tableA = tbls['YtblA']
# #tableB = tbls['YtblB']
# tableA = tbls['ZtblA']
# tableB = tbls['ZtblB']
# 
# plt.plot(tableA[0])
# plt.plot(tableA[-1])
# plt.show()
