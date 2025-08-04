from functools import partial
import tensorflow as tf
import random as rn
from scipy import optimize
import numpy as np
import sys
import time
import matplotlib.pyplot as plt
import scipy.io.matlab
import matplotlib
from DummyModel import DummyModel

from enum import Enum,unique
import keras
from keras import backend as K
from ANNmodel import AdvAutoencoder,datasetLoadUtility

#%% Options
class systemSelectorEnum(Enum):
    def loadFromDataset(filename,nonLinearInputChar=False):
        dynamicModel=DummyModel();
        dsLoading=datasetLoadUtility();
        Uvero, Yvero, UV,YV=dsLoading.loadDatasetFromMATfile(filename)        
        numel=Uvero.shape[0];
        numelV=UV.shape[0];
        u_n=np.reshape(Uvero.T[0],(numel,1))
        y_n=np.reshape(Yvero.T[0],(numel,1))
        u_Vn=np.reshape(UV.T[0],(numelV,1))
        y_Vn=np.reshape(YV.T[0],(numelV,1))
        #meanY=np.mean(y_n)
        #meanU=np.mean(u_n)
        #stdY=np.std(y_n)
        #stdU=np.std(u_n)
        #y_n=(y_n-meanY)/stdY;#+np.random.normal(0,0.05,(sizeT,1))
        #y_Vn=(y_Vn-meanY)/stdY;#+np.random.normal(0,0.05,(sizeV,1))
        #u_n=(u_n-meanU)/stdU;#+np.random.normal(0,0.05,(sizeT,1))
        #u_Vn=(u_Vn-meanU)/stdU;#+np.random.normal(0,0.05,(sizeV,1))
        return dynamicModel,u_n,y_n,u_Vn,y_Vn
    def wienerHammersteindataset():
        return systemSelectorEnum.loadFromDataset('datasets/wh.mat');
    def gyroscopedataset():
        return systemSelectorEnum.loadFromDataset('datasets/gyroscope.mat');
    def SILVERBOXdataset():
        return systemSelectorEnum.loadFromDataset('datasets/Silverbox.mat');


class Options():
    def __init__(self):    
        self.nonLinearInputChar=True; 
        self.dynamicalSystemSelector=systemSelectorEnum.wienerHammersteindataset
        self.stringDynamicalSystemSelector=str(self.dynamicalSystemSelector).replace('<function systemSelectorEnum.','').split(' at ')[0]
        self.affineStruct=True;
        self.openLoopStartingPoint=15
        self.horizon=5# for MPC
        self.TRsteps=1
        self.fitHorizon=5# 5 or 2 #it's +1 wrt to paper
        self.n_a=10;#n_a=n_b
        self.useGroupLasso=False;             
        self.stateReduction=True;
        self.regularizerWeight=0.0001;
        self.closedLoopSim=True
        self.enablePlot=False
        self.stateSize=6;
        pass

#%% Functions definition
def prepareMatrices(uSequence,x0):
    logY=[]
    logX=[]
    uSequence=np.array(uSequence)
    
    for u in uSequence:
        u=np.reshape(u,(1,1))
        x0=model.bridgeNetwork.predict([u,x0])   
        y=model.outputEncoder.predict([x0[0]])        
        logY+=[y]
        logX+=[x0]
        x0=x0[0]
        pass
    
    return logX,logY
    pass


def costFunction(uSequence,r,um1,logAB,logC,x0):
    
    logY=[]
    uSequence=np.array(uSequence)
    um1=np.array(um1)
    i=0
    for u in uSequence:
        #u=np.reshape(u,(1,1))
        asda=np.concatenate([x0.ravel(),u.ravel()])
        asda=np.reshape(asda,(Option.stateSize+1,1))
        x0=np.dot(logAB[i][1],asda)
        x0=x0+np.reshape(logAB[i][2],(Option.stateSize,1))
        y=np.dot(logC[i][1],x0)
        logY+=[y[0][-1]]
        i=i+1
        pass
#    logY+=[y[0][1]]
    logY=np.array(logY)
#    print(logY-r)
    cost=.001*np.sum(np.square(uSequence))+\
        .01*np.sum(np.square(uSequence[1:]-uSequence[:-1]))+\
        .01*np.sum(np.square(uSequence[0]-um1))+\
        np.sum(np.square(logY-r))*1
    
    return cost


def evaluateFeatureImportance():
    from matplotlib.ticker import MaxNLocator


    if not Option.stateReduction:
        w=model.convEncoder.get_layer('enc00').get_weights()
        ax = plt.figure(figsize=[8,2]).gca()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        neuronsCount=np.sum(abs(w[0])>1e-3,1);
#        print(len(neuronsCount))
        windowsLen=int(len(neuronsCount)/2)
        yAxis=range(0,windowsLen)[::-1]
        print(neuronsCount,'encoder=>')    
        plt.title('$encoder$')
        plt.step(yAxis,neuronsCount[0:windowsLen],where='mid')
        plt.step(yAxis,neuronsCount[windowsLen:],where='mid')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.tight_layout() 
        
    else:
        w1=model.bridgeNetwork.get_layer('bridge00').get_weights()
        w=model.outputEncoder.get_layer('dec00').get_weights()
        neuronsCount=np.sum(abs(w1[0][0:-1])>1e-3,1);
        yAxis=range(0,len(neuronsCount))
        print(neuronsCount,'bridge=>')  
        plt.figure(figsize=[8,2])
        plt.title('$bridge$')
        plt.step(yAxis,neuronsCount,where='mid')
        plt.tight_layout() 
        neuronsCount=np.sum(abs(w[0])>1e-3,1);        
        print(neuronsCount,'decoder=>')    
        yAxis=range(0,len(neuronsCount))
        ax = plt.figure(figsize=[8,2]).gca()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.title('$decoder$')
        plt.step(yAxis,neuronsCount,where='mid')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.tight_layout() 
    pass


def openLoopValidation(model,validationOnMultiHarmonic=True, reset=-1,YTrue=None,U_Vn=None,openLoopStartingPoint=15):    
    pastY=np.zeros((model.strideLen,1))
    pastU=np.zeros((model.strideLen,1))

    x0=model.convEncoder.predict([pastY.T,pastU.T])
    logY=[]
    logU=[]
    logYR=[]
    finalRange=1000;
    if not(YTrue is None):
        finalRange=YTrue.shape[0]
    for i in range(0,finalRange):
        y_kReal=YTrue[i]
        u=[U_Vn[i]]
        
        pastU=np.reshape(np.append(pastU,u)[1:],(model.strideLen,1))
        pastY=np.reshape(np.append(pastY,y_kReal)[1:],(model.strideLen,1))
        if i<openLoopStartingPoint or (i%reset==0 and reset>0):
            x0=model.convEncoder.predict([pastY.T,pastU.T])
            print('*',end='')
        else:        
            x0=model.bridgeNetwork.predict([u,x0])[0]
        y=model.outputEncoder.predict([x0])[0]
        if i>=openLoopStartingPoint:
            logY+=[y[0][-2]]
            logYR+=[y_kReal[0]]
            logU+=[u[0]]
        
        pass
    print('\n')
    logY=np.array(logY)
    logYR=np.array(logYR)    
    a=np.linalg.norm(np.array(logY)-np.array(logYR))
    b=np.linalg.norm(np.mean(np.array(logYR))-np.array(logYR))
    fit=1-(a/b)
    NRMSE=1-np.sqrt(np.mean(np.square(np.array(logY)-np.array(logYR))))/(np.max(logYR)-np.min(logYR))
    fit=np.max([0,fit])
    NRMSE=np.max([0,NRMSE])
    print('fit: ',fit)
    print('NRMSE: ',NRMSE)
    return fit,NRMSE,logY,logYR
