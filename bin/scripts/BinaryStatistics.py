#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from numpy import loadtxt
from scipy.stats import t

from GetShpFileInfo import getFieldElement
from GetShpFileInfo import getFieldElement,getCropCoverDictionary
from dates import GetDates

#--------------------------------------------------------------------

# TP (True Positive) : Crop Hit
# TN (True Negative) : Correct No-Crop rejection
# FP (False Positive) : It is a false alarm : No-Crop detected as crops
# FN (False Negative) : It is a miss : Crop detected as no-crops

def getBinaryStatistics(MatrixFileBase,Dates):  
    PreCrop = np.zeros(len(Dates))
    PreNoCrop = np.zeros(len(Dates))
    RecCrop = np.zeros(len(Dates))
    RecNoCrop = np.zeros(len(Dates))
    OA = np.zeros(len(Dates))
    for idDate,iDate in enumerate(Dates):
        MatrixFile =  MatrixFileBase + "_Date%i" % iDate +".csv"
        print MatrixFile
        ConfMat = loadtxt(MatrixFile, delimiter = ",")
        TN = ConfMat[0][0]; FN = ConfMat[0][1]
        FP = ConfMat[1][0];TP = ConfMat[1][1]
        PreCrop[idDate] = float(TP) / (FP + TP)
        RecCrop[idDate] = float(TP) / (FN + TP)
        PreNoCrop[idDate] = float(TN) / (FN + TN)
        RecNoCrop[idDate] = float(TN) / (FP + TN)
        OA[idDate] = float((TN + TP)) / (FP + TP + TN + FN)   
    return  PreCrop , RecCrop ,  PreNoCrop, RecNoCrop, OA


def getBinaryStatisticsSmart(MatrixFileBase,DatesMax):  
    PreCrop = np.zeros(len(DatesMax))
    PreNoCrop = np.zeros(len(DatesMax))
    RecCrop = np.zeros(len(DatesMax))
    RecNoCrop = np.zeros(len(DatesMAx))
    OA = np.zeros(len(DatesMax))
    for idDate,iDate in enumerate(DatesMax):
        MatrixFile =  MatrixFileBase + "_Date%i" % iDate +".csv"
        print MatrixFile
        ConfMat = loadtxt(MatrixFile, delimiter = ",")
        TN = ConfMat[0][0]; FN = ConfMat[0][1]
        FP = ConfMat[1][0];TP = ConfMat[1][1]
        PreCrop[idDate] = float(TP) / (FP + TP)
        RecCrop[idDate] = float(TP) / (FN + TP)
        PreNoCrop[idDate] = float(TN) / (FN + TN)
        RecNoCrop[idDate] = float(TN) / (FP + TN)
        OA[idDate] = float((TN + TP)) / (FP + TP + TN + FN)   
    return  PreCrop , RecCrop ,  PreNoCrop, RecNoCrop, OA



def getPerClassStatistics(CropValFile,CropPredFile,CropTypeValFile,CropTypePredFile):
    binval = loadtxt(CropValFile)
    binpred = loadtxt(CropPredFile)
    typeval = loadtxt(CropTypeValFile)
    typepred = loadtxt(CropTypePredFile)
    binpredfromtype = np.zeros(typepred.shape)
    
    classes,nbr,cropclasses,cropnbr = classcount(binval,typeval)
    for i in range(len(typepred)):
        if typepred[i] in cropclasses:
            binpredfromtype[i] = 1
        else:
            binpredfromtype[i] = 0
    
    classespred,nbrpred,cropclassespred,cropnbrpred = classcount(binpredfromtype,typepred)
    #classespred,nbrpred,cropclassespred,cropnbrpred = classcount(binpred,typepred) Does not work because some classe appears sooner 
    #                                                                               in the prediction file wich create an error.
    #                                                                               Good way above even if it might look undirect and useless.

    stat = []
    # 1. PC = Class Precision    
    # 2. PosRC = Positive recall
    # 3. NegRC = Negative recall
    # 4. LNegP = Lost Negative Precision
    # 5. LPosP = Lost Positive Precision
    for i in range(len(classes)):
        PC = 0
        PosRC = 0
        NegRC = 0
        LNegP = 0
        LPosP = 0
        ci = int(classes[i])			# ci = cp
        Ncr = nbr[classes.index(ci)]		# Nbr of ci element in validation set 
        Ncp = nbrpred[classespred.index(ci)]	# Nbr of ci element in prediction set 
        for j in range(len(binval)):		# Loop on validation set
            r = int(binval[j])
            cr = int(typeval[j])
            cp = int(typepred[j])
            if ci in cropclasses:pi = 1
            else:pi = 0
            if cr in cropclasses:pr = 1
            else:pr = 0
            if cp in cropclasses:pp = 1
            else:pp = 0

            if ci == cp:
                if cp == cr: PC += 1			# Checked
                if cp != cr and pp == pr: LPosP += 1	# Checked
                if cp != cr and pp != pr: LNegP += 1	# Checked
          
            if ci == cr:
                if cr != cp and pr == pp: PosRC += 1	# Checked
                if cr != cp and pr != pp: NegRC += 1    # Checked
        
        RC = PC     
        PC = 100*float(PC)/float(Ncp)
        RC = 100*float(RC)/float(Ncr)
        LNegP = 100*float(LNegP)/float(Ncp)
        LPosP = 100*float(LPosP)/float(Ncp)
        PosRC = 100*float(PosRC)/float(Ncr)
        NegRC = 100*float(NegRC)/float(Ncr)
        stat.append([PC,LPosP,LNegP,RC,PosRC,NegRC])  
    return stat

def getPerClassStatisticsInduced(CropValFile,CropTypeValFile,CropTypePredFile):
    binval = loadtxt(CropValFile)
    typeval = loadtxt(CropTypeValFile)
    typepred = loadtxt(CropTypePredFile)
    binpred = np.zeros(typepred.shape)
    
    classes,nbr,cropclasses,cropnbr = classcount(binval,typeval)
    for i in range(len(typepred)):
        if typepred[i] in cropclasses:
            binpred[i] = 1
        else:
            binpred[i] = 0
    
    classespred,nbrpred,cropclassespred,cropnbrpred = classcount(binpred,typepred)
    stat = []
    # N1 = [ 0 for i in range(len(classes))]
    # 1. PC = Class Precision    
    # 2. PosRC = Positive recall
    # 3. NegRC = Negative recall
    # 4. LNegP = Lost Negative Precision
    # 5. LPosP = Lost Positive Precision
    #    Ncrop = Total of crop (minus considered class if it is a crop)
    #    Nnocrop = Total of no crop (minus considered class if it is a nocrop)
    for i in range(len(classes)):
        PC = 0
        PosRC = 0
        NegRC = 0
        LNegP = 0
        LPosP = 0
        ci = int(classes[i])			# ci = cp
        Ncr = nbr[classes.index(ci)]		# Nbr of ci element in validation set 
        Ncp = nbrpred[classespred.index(ci)]	# Nbr of ci element in prediction set 
        for j in range(len(binval)):		# Loop on validation set
            r = int(binval[j])
            cr = int(typeval[j])
            cp = int(typepred[j])
            if ci in cropclasses:pi = 1
            else:pi = 0
            if cr in cropclasses:pr = 1
            else:pr = 0
            if cp in cropclasses:pp = 1
            else:pp = 0

            if ci == cp:
                if cp == cr: PC += 1			# Checked
                if cp != cr and pp == pr: LPosP += 1	# Checked
                if cp != cr and pp != pr: LNegP += 1	# Checked
          
            if ci == cr:
                if cr != cp and pr == pp: PosRC += 1	# Checked
                if cr != cp and pr != pp: NegRC += 1    # Checked
        RC = PC     
        PC = 100*float(PC)/float(Ncp)
        RC = 100*float(RC)/float(Ncr)
        LNegP = 100*float(LNegP)/float(Ncp)
        LPosP = 100*float(LPosP)/float(Ncp)
        PosRC = 100*float(PosRC)/float(Ncr)
        NegRC = 100*float(NegRC)/float(Ncr)
        stat.append([PC,LPosP,LNegP,RC,PosRC,NegRC])  
    return stat


def CI95(x):
   df = len(x) - 1
   q = t.ppf(0.95,df)
   sigma = np.std(x,axis = 0)
   CI = q*sigma/np.sqrt(df+1) 
   return CI

def classcount(y,x): # Attention x et y inversé
   classes = []
   nbr = []
   cropclasses = []
   cropnbr = []
   for i in range(len(x)):
      cl = int(x[i])
      iscrop = int(y[i])
      if cl not in classes:
         classes.append(cl)  
         nbr.append(1)
         if iscrop == 1:
             cropclasses.append(cl)  
             cropnbr.append(1)
      else:
          nbr[classes.index(cl)] += 1
          if iscrop == 1:
              cropnbr[cropclasses.index(cl)] += 1


   z = np.array((classes,nbr))
   z = z.transpose()
   z = np.array(sorted(z, key=lambda row: row[0]))
   classes = z.transpose()[0].tolist()
   nbr = z.transpose()[1].tolist()
   
   z = np.array((cropclasses,cropnbr))
   z = z.transpose()
   z = np.array(sorted(z, key=lambda row: row[0]))
   cropclasses = z.transpose()[0].tolist()
   cropnbr = z.transpose()[1].tolist()
   
   return classes,nbr,cropclasses,cropnbr

#--------------------------------------------------------------------

btype = "NDWI_"
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
    BinaryFile = DirRun + "CropLabels_val.txt"
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
    print "* Direct Binary Classification for Run_%i *" % (tirage) 
    print "Number of Crop pixels:\t\t",NbrCrop
    print "Number of NoCrop pixels:\t",NbrNoCrop
    print "Total:\t\t\t\t",NbrCrop + NbrNoCrop
  
    Dates = []
    for iDate in range(1, NbDates+1):
        Dates.append (iDate)
    PreCrop , RecCrop ,  PreNoCrop, RecNoCrop, OverallA = getBinaryStatistics(DirRun+"RF_ResultsBinary/ConfuMatri/" + btype + "ConfuMatrix",Dates)
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
  print ""
  print "Precision and recall per class"

  AverageStatPerClass = [] 
  for tirage in range(Nbtirages):
    print "Run ",tirage
    DirRun = InDir + "Run_%i/" % (tirage)
    cropvalfile = DirRun + "/CropLabels_val.txt"
    croppredfile = DirRun + "RF_ResultsBinary/Predicted/NDWI_PredictedLabels_Date%d.txt"%(NbDates)
    typevalfile = DirRun + "/CropTypeLabels_val.txt"
    typepredfile = DirRun + "RF_ResultsAll/Predicted/NDWI_PredictedLabels_Date%d.txt"%(NbDates)

    statperclass = getPerClassStatistics(cropvalfile,croppredfile,typevalfile,typepredfile)
    AverageStatPerClass.append(statperclass)

  mat = np.asarray(AverageStatPerClass)
  MoyStat = np.mean(mat,axis=0)
  print MoyStat

