import matplotlib.pyplot as plt
import numpy as np

from smt.applications.mixed_integer import MixedIntegerKrigingModel

from smt.surrogate_models import KRG
from smt.design_space import (
    DesignSpace,
    IntegerVariable,
)


sm_q = KRG.load("SM_heatFlux")


fileName  = '../YValidation.dat'
fileName2 = '../Validation_Xdata.dat'


XV = []
YV = []

with open(fileName) as file:
    file2 = open(fileName2)
    lines  = file.readlines()
    lines2 = file2.readlines()
    
    line_index = 0
    for line in lines:
        if line_index > 0:
            splitted = line.split()
            if splitted[2] == 'NotConv':
                splitted[2] = 999
            conv = float(splitted[2])
            if conv < 1e-2:
                YV.append(float(splitted[1]))
                x = lines2[line_index].split()
                XV.append([float(x[1]), float(x[2]), float(x[3]), float(x[4]), float(x[5])])
        line_index += 1
file2.close()

XV = np.array(XV)
YV = np.array(YV)

YV_K = sm_q.predict_values(XV)

plt.scatter(YV[:], YV_K[:])
plt.plot([np.min(YV[:]), np.max(YV[:])], [np.min(YV[:]), np.max(YV[:])])
plt.show()
