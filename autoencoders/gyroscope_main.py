import numpy as np
import time
import matplotlib.pyplot as plt

np.random.seed(1)
from ANNmodel import AdvAutoencoder
from functions import *

Option = Options()
Option.fitHorizon = 11
Option.dynamicalSystemSelector = systemSelectorEnum.gyroscopedataset
Option.closedLoopSim = False
Option.stringDynamicalSystemSelector = str( \
    Option.dynamicalSystemSelector).replace('<function systemSelectorEnum.',
    '').split(' at ')[0]
Option.nonLinearInputChar = False
Option.stateSize = 5
Option.n_a = 10
Option.affineStruct = True
Option.useGroupLasso = False;   
Option.regularizerWeight = 0.0001

import warnings
warnings.filterwarnings("ignore")

#%% DS generation and model learning
simulatedSystem, U_n, Y_n, U_Vn, Y_Vn = Option.dynamicalSystemSelector()

model = AdvAutoencoder(affineStruct = Option.affineStruct,
        useGroupLasso = Option.useGroupLasso,
        stateReduction = Option.stateReduction,
        fitHorizon = Option.fitHorizon,
        strideLen = Option.n_a,#n_a=n_b
        outputWindowLen = 2,#+1 wrt the paper
        n_layer = 2,
        n_neurons = 64,                     
        regularizerWeight = Option.regularizerWeight,
        stateSize = Option.stateSize)
model.setDataset(U_n.copy(), Y_n.copy(), U_Vn.copy(), Y_Vn.copy())

inputU, inputY = model.prepareDataset()
model.trainModel()
predictedLeft, stateLeft, oneStepAheadPredictionError, \
    forwardedPredictedError, forwardError = \
    model.model.predict([inputY, inputU])

#%% Model Validation Validation
voM = False
r = -1
fit, NRMSE, logY, logYR = openLoopValidation(model,
        validationOnMultiHarmonic = voM,
        reset = r,
        YTrue = Y_Vn.copy(),
        U_Vn = U_Vn.copy(),
        openLoopStartingPoint = Option.openLoopStartingPoint)

r = 10
fit, NRMSE, logY, logYR = openLoopValidation(model,
        validationOnMultiHarmonic = voM,
        reset = r,
        YTrue = Y_Vn.copy(),
        U_Vn = U_Vn.copy(),
        openLoopStartingPoint = Option.openLoopStartingPoint)

k = 1000
plt.figure()
plt.plot(logYR[1:k], label = 'y')
plt.plot(logY[1:k], label = 'hy')
plt.legend()
plt.show()
