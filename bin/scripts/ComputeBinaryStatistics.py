#!/usr/bin/python
# -*- coding: utf-8 -*-
### CesBIO 2019 ###
import os
import sys
import numpy as np
from numpy import loadtxt
#from scipy.stats import t

from GetShpFileInfo import getFieldElement
from GetShpFileInfo import getFieldElement,getCropCoverDictionary
from dates import GetDates

#--------------------------------------------------------------------

# TP (True Positive) : Crop Hit
# TN (True Negative) : Correct No-Crop rejection
# FP (False Positive) : It is a false alarm : No-Crop detected as crops
# FN (False Negative) : It is a miss : Crop detected as no-crops
def getDirectBinaryStatistics(MatrixFileBase,Dates):  
    PreCrop = np.zeros(len(Dates))
    PreNoCrop = np.zeros(len(Dates))
    RecCrop = np.zeros(len(Dates))
    RecNoCrop = np.zeros(len(Dates))
    OA = np.zeros(len(Dates))
    for idDate,iDate in enumerate(Dates):
        MatrixFile =  MatrixFileBase + "_Date%i" % iDate +".csv"
        print MatrixFile
        ConfMat = loadtxt(MatrixFile, delimiter = ",")
        TN = ConfMat[0][0]; FP = ConfMat[0][1]
        FN = ConfMat[1][0];TP = ConfMat[1][1]
        PreCrop[idDate] = float(TP) / (FP + TP)
        RecCrop[idDate] = float(TP) / (FN + TP)
        PreNoCrop[idDate] = float(TN) / (FN + TN)
        RecNoCrop[idDate] = float(TN) / (FP + TN)
        OA[idDate] = float((TN + TP)) / (FP + TP + TN + FN)   
    return  PreCrop , RecCrop ,  PreNoCrop, RecNoCrop, OA

def getBinaryStatistics(PredictedFile,ReferenceFile,ShapeFile,lc,code,crop):
    PredictedSet = loadtxt(PredictedFile, delimiter=',')
    ReferenceSet = loadtxt(ReferenceFile, delimiter=',')

    if len(PredictedSet) != len(ReferenceSet):
        raise Exception(" Predicted and reference data set do not contain the same number of samples")
    
    CropDict = getCropCoverDictionary(ShapeFile,lc,code,crop)

    TP = TN = FP = FN = 0
    for id,ReferenceSample in enumerate(ReferenceSet):
        # If the reference sample is a crop
        if PredictedSet[id] == 0:
          print id,ReferenceSample,PredictedSet[id]
          quit()
        if CropDict[ReferenceSample] == 1:
            #print CropDict[PredictedSet[id]]
            if CropDict[ReferenceSample] == CropDict[PredictedSet[id]]:
                TP = TP + 1
            else:
                FN = FN + 1
        
        # If the reference sample is a non crop
        else:
            if CropDict[ReferenceSample] == CropDict[PredictedSet[id]]:
                TN = TN + 1
            else:
                FP = FP + 1

    return TP, TN, FP, FN

#--------------------------------------------------------------------


def getBinaryQualityMeasures(PredictedFile,ReferenceFile,ShapeFile,lc,code,crop):

    [TP,TN,FP,FN] = getBinaryStatistics(PredictedFile,ReferenceFile,ShapeFile,lc,code,crop)
    
    PreCrop = float(TP) / (FP + TP)

    RecCrop = float(TP) / (FN + TP)

    PreNoCrop = float(TN) / (FN + TN)

    RecNoCrop = float(TN) / (FP + TN)

    OA = float((TN + TP)) / (FP + TP + TN + FN)

    return  PreCrop , RecCrop ,  PreNoCrop, RecNoCrop, OA

#--------------------------------------------------------------------

def getBinaryTemporalStatistics(InPredictedDir,btype,ReferenceFile,Dates,ShapeFile,lc,code,crop):
    
    OA = np.zeros(len(Dates))
    PreCrop = np.zeros(len(Dates))
    PreNoCrop = np.zeros(len(Dates))
    RecCrop = np.zeros(len(Dates))
    RecNoCrop = np.zeros(len(Dates))

    for idDate,iDate in enumerate(Dates):
        print "date =",idDate 
        PredictedFile =  InPredictedDir + btype +"PredictedLabels_Date%i" % iDate +".txt"
        [PreCrop[idDate] , RecCrop[idDate] , PreNoCrop[idDate], RecNoCrop[idDate], OA[idDate]] = getBinaryQualityMeasures(PredictedFile,ReferenceFile,ShapeFile,lc,code,crop)

    return  PreCrop , RecCrop ,  PreNoCrop, RecNoCrop, OA

def CI95(x):
        tdistOneSided={4:1.132,9:1.833}
        tdistTwoSided={4:2.776,9:2.262}
        df = len(x) - 1
        #q = tdist.ppf(0.95,df)
        q = tdistOneSided[df]
        sigma = np.std(x,axis = 0)
        CI = q*sigma/np.sqrt(df+1) 
        return CI

