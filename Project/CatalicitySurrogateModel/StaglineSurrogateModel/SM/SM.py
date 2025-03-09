import matplotlib.pyplot as plt
import numpy as np

from smt.applications.mixed_integer import MixedIntegerKrigingModel

from smt.surrogate_models import KRG
from smt.design_space import (
    DesignSpace,
    IntegerVariable,
)

fileName  = '../YTraining.dat'
fileName2 = '../Training_Xdata.dat'

XT = []
YT = []
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
                YT.append(float(splitted[1]))
                x = lines2[line_index].split()
                XT.append([float(x[1]), float(x[2]), float(x[3]), float(x[4]), float(x[5])])
        line_index += 1
file2.close()
print(f'{XT[0]} {YT[0]}')
        
fileName  = '../YValidation.dat'
fileName2 = '../Validation_Xdata.dat'

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

XT = np.array(XT)
YT = np.array(YT)
XV = np.array(XV)
YV = np.array(YV)

design_space = DesignSpace(
    [
        IntegerVariable(0, 1),
    ]
)

sm_q = KRG(theta0=[1e-1], corr='matern52')
#sm_drag.corr = 'matern52'
#sm_drag = MixedIntegerKrigingModel(surrogate=KRG(design_space=design_space, theta0=[1e-2], hyper_opt="Cobyla"))
sm_q.set_training_values(XT, YT)
sm_q.train()
sm_q.save('./SM_heatFlux')

YV_K = sm_q.predict_values(XV)
err = 0
for i in range(len(YV_K)):
    err += (YV[i] - YV_K[i])**2
    #print(f'{YV[i,0]} {YV_K[i]}')
err = np.sqrt(err)/i*100/np.mean(YV[i])
print(err)

plt.scatter(YV[:], YV_K[:])
plt.plot([np.min(YV[:]), np.max(YV[:])], [np.min(YV[:]), np.max(YV[:])])
plt.show()
