#!/usr/bin/python
# -*- coding: utf-8 -*-
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

def getBinaryStatistics(PredictedFile,ReferenceFile,ShapeFile):
    PredictedSet = loadtxt(PredictedFile, delimiter=',')
    ReferenceSet = loadtxt(ReferenceFile, delimiter=',')

    if len(PredictedSet) != len(ReferenceSet):
        raise Exception(" Predicted and reference data set do not contain the same number of samples")
    
    CropDict = getCropCoverDictionary(ShapeFile)

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


def getBinaryQualityMeasures(PredictedFile,ReferenceFile,ShapeFile):

    [TP,TN,FP,FN] = getBinaryStatistics(PredictedFile,ReferenceFile,ShapeFile)
    
    PreCrop = float(TP) / (FP + TP)

    RecCrop = float(TP) / (FN + TP)

    PreNoCrop = float(TN) / (FN + TN)

    RecNoCrop = float(TN) / (FP + TN)

    OA = float((TN + TP)) / (FP + TP + TN + FN)

    return  PreCrop , RecCrop ,  PreNoCrop, RecNoCrop, OA

#--------------------------------------------------------------------


def getBinaryTemporalStatistics(InPredictedDir,btype,ReferenceFile,Dates,ShapeFile):
    
    OA = np.zeros(len(Dates))
    PreCrop = np.zeros(len(Dates))
    PreNoCrop = np.zeros(len(Dates))
    RecCrop = np.zeros(len(Dates))
    RecNoCrop = np.zeros(len(Dates))

    for idDate,iDate in enumerate(Dates):
        print "date =",idDate 
        PredictedFile =  InPredictedDir + btype +"PredictedLabels_Date%i" % iDate +".txt"
        [PreCrop[idDate] , RecCrop[idDate] , PreNoCrop[idDate], RecNoCrop[idDate], OA[idDate]] = getBinaryQualityMeasures(PredictedFile,ReferenceFile,ShapeFile)

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

#--------------------------------------------------------------------


if __name__ == "__main__":

  dates,datesLabels = GetDates("/home/arnaudl/Documents/OTB/Sentinel2-Data/S2_DateListF.txt")
  NbDates = 33
  Nbtirages = 5
  InDir = "OpticalSamples/"
  PC = np.zeros((Nbtirages,NbDates))
  RC = np.zeros((Nbtirages,NbDates))
  PNoC = np.zeros((Nbtirages,NbDates))
  RNoC = np.zeros((Nbtirages,NbDates))
  OA = np.zeros((Nbtirages,NbDates))
  for tirage in range(Nbtirages):
    DirRun = InDir + "Run_%i/" % (tirage) 
    PredictedDir = DirRun + "RF_ResultsAll/Predicted/"
    ReferenceFile = DirRun + "CropTypeLabels_val.txt"
    BinaryFile = DirRun + "CropLabels_val.txt"
    ShapeFile = "../../ShapeFiles/CombineFinal-L93-20Merosion_Ludovic.shp"

    ReferenceBinary = loadtxt(BinaryFile, delimiter=',')
    NbrCrop = 0
    NbrNoCrop = 0
    for value in ReferenceBinary:
      if int(value) == 0:
        NbrNoCrop += 1
      else:
        NbrCrop += 1
  
    print ""
    print "--------------------------------------------"
    print "* Deduced Binary Classification for Run_%i *" % (tirage) 
    print "Number of Crop pixels:\t\t",NbrCrop
    print "Number of NoCrop pixels:\t",NbrNoCrop
    print "Total:\t\t\t\t",NbrCrop + NbrNoCrop
  
    Dates = []
    for iDate in range(1, NbDates+1):
        Dates.append (iDate)
    print PredictedDir
    print ReferenceFile
    print BinaryFile

    PreCrop , RecCrop ,  PreNoCrop, RecNoCrop, OverallA = getBinaryTemporalStatistics(PredictedDir,ReferenceFile,Dates,ShapeFile)
    PC[tirage] = 100*PreCrop
    RC[tirage] = 100*RecCrop
    PNoC[tirage] = 100*PreNoCrop
    RNoC[tirage] = 100*RecNoCrop
    OA[tirage] = 100*OverallA
  aPC = np.mean(PC,axis = 0)
  aRC = np.mean(RC,axis = 0)
  aPNoC = np.mean(PNoC,axis = 0)
  aRNoC = np.mean(RNoC,axis = 0)
  aOA = np.mean(OA,axis = 0)
  sigma = np.std(OA,axis = 0)
    
  CI = CI95(OA)   
  print "Date\tPreC\tRecC\tPreNoC\tRecNoC\tOA\tCI95"
  print "-----------------------------------------------------"
  for d in range(len(PreCrop)):
    print "%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f" %(datesLabels[d],aPC[d],aRC[d],aPNoC[d],aRNoC[d],aOA[d],CI[d])
  print "-----------------------------------------------------"
