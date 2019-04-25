#!/usr/bin/env python
#-*- coding: utf-8 -*-

###################################################################################################
# Plot automatical the average fscore for optical, radar or both depending on the input parameter #
# Ludo 14/06/2017  STILL IN PROGRESS. SO FAR, IT WORKS WELL BUT IT IS REALLY NOT ELEGANT          #
###################################################################################################

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'tolatex'))
from collections import defaultdict
from scipy.stats import t as tdist
from numpy import loadtxt
from dates import GetDates,GetDatesSentinel
import numpy as np
import argparse
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import csv
import math
import tolatex as latex
import ComputeBinaryStatistics as BinStat
import BinaryStatistics as ClassBinStat


from osgeo import ogr, osr
try:
  from osgeo import gdal
except ImportError:
  import gdal


def CI95(x):
	tdistOneSided={1:6.314,4:2.132,8:1.860,9:1.833}
	tdistTwoSided={1:12.706,4:2.776,8:2.306,9:2.262}
        df = len(x) - 1
        #q = tdist.ppf(0.95,df)
        q = tdistOneSided[df]
        sigma = np.std(x,axis = 0)
        CI = q*sigma/np.sqrt(df+1) 
        return CI

def computeKappa(confMat):

	nbrGood = confMat.trace()
	nbrSample = confMat.sum()

	overallAccuracy  = float(nbrGood) / float(nbrSample)

	## the lucky rate.
	luckyRate = 0.
	for i in range(0, confMat.shape[0]):
		sum_ij = 0.
       		sum_ji = 0.
        	for j in range(0, confMat.shape[0]):
         		sum_ij += confMat[i][j]
                	sum_ji += confMat[j][i]
        	luckyRate += sum_ij * sum_ji

	# Kappa.
	if float((nbrSample*nbrSample)-luckyRate) != 0:
		kappa = float((overallAccuracy*nbrSample*nbrSample)-luckyRate)/float((nbrSample*nbrSample)-luckyRate)
	else :
		kappa = 1000

	return kappa

def computePreByClass(confMat,AllClass):

	Pre = []#[(class,Pre),(...),()...()]

	for i in range(len(AllClass)):
		denom = 0
		for j in range(len(AllClass)):
			denom += confMat[j][i]
			if i == j:
				nom = confMat[j][i]
		if denom != 0:
			currentPre = float(nom)/float(denom)
		else :
			currentPre = 0.
		Pre.append((AllClass[i],currentPre))
	return Pre

def computeRecByClass(confMat,AllClass):
	Rec = []#[(class,rec),(...),()...()]
	for i in range(len(AllClass)):
		denom = 0
		for j in range(len(AllClass)):
			denom += confMat[i][j]
			if i == j:
				nom = confMat[i][j]
		if denom != 0 :
			currentRec = float(nom)/float(denom)
		else:
			currentRec = 0.
		Rec.append((AllClass[i],currentRec))
	return Rec

def computeFsByClass(Pre,Rec,AllClass):
	Fs = []
	for i in range(len(AllClass)):
		if float(Rec[i][1]+Pre[i][1]) != 0:
			Fs.append((AllClass[i],float(2*Rec[i][1]*Pre[i][1])/float(Rec[i][1]+Pre[i][1])))
		else:
			Fs.append((AllClass[i],0.0))
	return Fs

def confCoordinatesCSV(csvPath):
	"""
	IN :
		csvPaths [string] : path to csv file
			ex : "/path/to/file1.csv"
	OUT : 
		out [list of lists] : containing csv's coordinates

		ex : file1.csv
			#Reference labels (rows):11
			#Produced labels (columns):11,12
			14258,52

		     file2.csv
			#Reference labels (rows):12
			#Produced labels (columns):11,12
			38,9372

		out = [[12,[11,38]],[12,[12,9372]],[11,[11,14258]],[11,[12,52]]]
	"""
	out = []
	cpty = 0
	FileMat = open(csvPath,"r")
	while 1:
		data = FileMat.readline().rstrip('\n\r')
		if data == "":
			FileMat.close()
			break
		if data.count('#Reference labels (rows):')!=0:
			ref = data.split(":")[-1].split(",")
		elif data.count('#Produced labels (columns):')!=0:
			prod = data.split(":")[-1].split(",")
		else:
			y = ref[cpty]
			line = data.split(",")
			cptx = 0
			for val in line:
				x = prod[cptx]
				out.append([int(y),[int(x),float(val)]])
				cptx+=1
			cpty +=1
	return out

def gen_confusionMatrix(csv_f,AllClass,AllClass_prod):

	NbClasses = len(AllClass)

	confMat = [[0]*NbClasses]*NbClasses
	confMat = np.asarray(confMat)
	
	row = 0
	for classRef in AllClass:
		flag = 0#in order to manage the case "this reference label was never classified"
		for classRef_csv in csv_f:
			if classRef_csv[0] == classRef:
				col = 0
				for classProd in AllClass:
					for classProd_csv in classRef_csv[1]:
						if classProd_csv[0] == classProd:
							confMat[row][col] = confMat[row][col] + classProd_csv[1]
					col+=1
		row+=1
	return confMat



def getAllClass(cMatrix):

	with open(cMatrix) as f:
		for line in f:
        		if "#Reference labels" in line:
				AllClass_ref = line.rstrip('\n\r').split(":")[-1].split(",")
				AllClass_ref = [int(Currentclass) for Currentclass in AllClass_ref]
				nbRef = len(AllClass_ref)
			elif "#Produced labels (columns)" in line:
				AllClass_prod = line.rstrip('\n\r').split(":")[-1].split(",")
				AllClass_prod = [int(Currentclass) for Currentclass in AllClass_prod]
				nbProd = len(AllClass_prod)
		#if nbRef!=nbProd:
		#	print AllClass_prod
		#	print AllClass_ref
			#for classRef in AllClass_ref:
			 #     for classProd in AllClass_prod:
			  #	 if classRef == classProd:
			#		print classRef
			#raise Exception(cMatrix+" is not a square matrix")
		return AllClass_ref,AllClass_prod
		




class QualityMesureStruct(object):
	def __init__(self, OA, Kappa, Pre, Rec, Fs):
		self.OA = OA
		self.Kappa = Kappa
		self.Pre = Pre
		self.Rec =  Rec
		self.Fs =  Fs

#cMatrix is the path of the confusion matrix
def genCoeff(cMatrix,normalized=False):
	[AllClass,AllClass_prod] = getAllClass(cMatrix)

	mat = confCoordinatesCSV(cMatrix)
	d = defaultdict(list)

	for k,v in mat:
		d[k].append(v)

	csv_f = list(d.items())

        if(normalized):
		mat_sort_nn = gen_confusionMatrix(csv_f,AllClass,AllClass_prod)
		mat_sort = mat_sort_nn.astype(np.float)
                for i,row in enumerate(mat_sort):
			  mat_sort[i] = row/np.sum(row)
                #print mat_sort_nn
		#print mat_sort
	        #print
	else:
		mat_sort = gen_confusionMatrix(csv_f,AllClass,AllClass_prod)

        nbrGood = mat_sort.trace()
	nbrSample = mat_sort.sum()
	overallAccuracy  = float(nbrGood) / float(nbrSample)


	Kappa = computeKappa(mat_sort)
	Pre = computePreByClass(mat_sort,AllClass)
	Rec = computeRecByClass(mat_sort,AllClass)
	FS = computeFsByClass(Pre,Rec,AllClass)

	return QualityMesureStruct(overallAccuracy,Kappa, Pre, Rec, FS)

def genCoeffDirect(mat_sort,AllClass):
	nbrGood = mat_sort.trace()
	nbrSample = mat_sort.sum()
	overallAccuracy  = float(nbrGood) / float(nbrSample)
	Kappa = computeKappa(mat_sort)
	Pre = computePreByClass(mat_sort,AllClass)
	Rec = computeRecByClass(mat_sort,AllClass)
	FS = computeFsByClass(Pre,Rec,AllClass)
	return QualityMesureStruct(overallAccuracy,Kappa, Pre, Rec, FS)



def genCoeffGather(cMatrix,gatheringSet):

        # Create path for CROP confusion matrix from ALL confusion matrix path
        #pos = cMatrix.find("ConfusionsMatrix_") + len("ConfusionsMatrix_")
        #cMatrixCrop = cMatrix[:pos] + "CROP_" + cMatrix[pos:]
        # Get classes list for both cases
        #[AllClass,AllClass_prod] = getAllClass(cMatrix)

        [CropClass,CropClass_prod] = getAllClass(cMatrix)
 
        # construct new class sets 
        #cropset = set(CropClass)
        #ClassesToGather = [c for c in AllClass if c not in cropset]
        #NewClass = [0] + CropClass

        # Get ALL and CROP confusion Matrix
	#mat = confCoordinatesCSV(cMatrix)
	#d = defaultdict(list)
	#for k,v in mat:
        #	d[k].append(v)

	#csv_f = list(d.items())
	#mat_sort = gen_confusionMatrix(csv_f,AllClass,AllClass_prod)

	matCrop = confCoordinatesCSV(cMatrix)
        dc = defaultdict(list)
	for k,v in matCrop:
		dc[k].append(v)

	csv_f_crop = list(dc.items())
	#mat_sort = gen_confusionMatrix(csv_f,AllClass,AllClass_prod)
	mat_sort_crop = gen_confusionMatrix(csv_f_crop,CropClass,CropClass_prod)

        # Creat GatherIdx
        maxlabel = max(gatheringSet.itervalues())

        GatherIdx = [[] for i in range(maxlabel)]
        for i,c in enumerate(CropClass):
          GatherIdx[gatheringSet[c]-1].append(i) 

        # Calculate new confusion matrix
        gatherMat = CalculateMatrixMultiple(mat_sort_crop,GatherIdx)
    
	NewClass = [(i+1) for i in range(len(GatherIdx))]

        return genCoeffDirect(gatherMat,NewClass),gatherMat



def genCoeffImperfect(cMatrix):
        # Create path for CROP confusion matrix from ALL confusion matrix path
        pos = cMatrix.find("ConfusionsMatrix_") + len("ConfusionsMatrix_")
        cMatrixCrop = cMatrix[:pos] + "CROP_" + cMatrix[pos:]

        # Get classes list for both cases
        [AllClass,AllClass_prod] = getAllClass(cMatrix)
        [CropClass,CropClass_prod] = getAllClass(cMatrixCrop)

        # construct new class sets 
        cropset = set(CropClass)
        ClassesToGather = [c for c in AllClass if c not in cropset]
        NewClass = [0] + CropClass

        # Get ALL and CROP confusion Matrix
	mat = confCoordinatesCSV(cMatrix)
	matCrop = confCoordinatesCSV(cMatrixCrop)
	d = defaultdict(list)
	for k,v in mat:
		d[k].append(v)

	csv_f = list(d.items())
	mat_sort = gen_confusionMatrix(csv_f,AllClass,AllClass_prod)

        dc = defaultdict(list)
	for k,v in matCrop:
		dc[k].append(v)

	csv_f_crop = list(dc.items())
	mat_sort = gen_confusionMatrix(csv_f,AllClass,AllClass_prod)
	mat_sort_crop = gen_confusionMatrix(csv_f_crop,CropClass,CropClass_prod)

        # Found index for gathering
        idx1 = []
        idx2 = []
        pos = 0
#        #print "pos=",pos
        for i,c in enumerate(AllClass):
            if(c in ClassesToGather):
                idx1.append(i) 
            else:
                idx2.append(i)

        # Calculate new confusion matrix
        gatherMat = CalculateMatrix(mat_sort,idx1,idx2)
  
        ## Method 1
        newmat = np.copy(gatherMat)
        for i in range(1,len(gatherMat)):
          for j in range(1,len(gatherMat[0])):
            newmat[i,j] = mat_sort_crop[i-1,j-1]

        #print "All:"
        #print mat_sort
        if(False):
            print "CT perfect CM:"
            print mat_sort_crop
            print "All with No Crop gathering:"
            print gatherMat
            print "CP + NoCrop:"
            print newmat
            print "********************"

#       # Place reshuffle hear if required later. Not need so far     
#        ori_idx = range(len(idx2)+1)
#        new_idx = range(1,pos+1)+[0]+range(pos+1,len(idx2)+1)
#        #print "ori_idx = ",ori_idx
#        #print "new_idx = ",new_idx
#    
#        newmat[:,ori_idx] = newmat[:,new_idx]
#        newmat[ori_idx,:] = newmat[new_idx,:]
#    
        return genCoeffDirect(newmat,NewClass),newmat

def genCoeffDirectCT(cMatrix):
        # Create path for CROP confusion matrix from ALL confusion matrix path
        pos = cMatrix.find("ConfusionsMatrix_") + len("ConfusionsMatrix_")
        cMatrixCrop = cMatrix[:pos] + "CROP_" + cMatrix[pos:]

        # Get classes list for both cases
        [AllClass,AllClass_prod] = getAllClass(cMatrix)
        [CropClass,CropClass_prod] = getAllClass(cMatrixCrop)

        # construct new class sets 
        cropset = set(CropClass)
        ClassesToGather = [c for c in AllClass if c not in cropset]
        NewClass = [0] + CropClass

        # Get ALL and CROP confusion Matrix
	mat = confCoordinatesCSV(cMatrix)
	matCrop = confCoordinatesCSV(cMatrixCrop)
	d = defaultdict(list)
	for k,v in mat:
		d[k].append(v)

	csv_f = list(d.items())
	mat_sort = gen_confusionMatrix(csv_f,AllClass,AllClass_prod)

        dc = defaultdict(list)
	for k,v in matCrop:
		dc[k].append(v)

	csv_f_crop = list(dc.items())
	mat_sort = gen_confusionMatrix(csv_f,AllClass,AllClass_prod)
	mat_sort_crop = gen_confusionMatrix(csv_f_crop,CropClass,CropClass_prod)

        # Found index for gathering
        idx1 = []
        idx2 = []
        pos = 0
#        #print "pos=",pos
        for i,c in enumerate(AllClass):
            if(c in ClassesToGather):
                idx1.append(i) 
            else:
                idx2.append(i)

        # Calculate new confusion matrix
        gatherMat = CalculateMatrix(mat_sort,idx1,idx2)
  
        # Method 2
        newmat = np.copy(mat_sort_crop)
        for i in range(len(mat_sort_crop)):
          for j in range(len(mat_sort_crop[0])):
            newmat[i,j] = gatherMat[i+1,j+1]


        #print "All:"
        #print mat_sort
        if(False):
            print "CT perfect CM:"
            print mat_sort_crop
            print "All with No Crop gathering:"
            print gatherMat
            print "CP + NoCrop:"
            print newmat
            print "********************"

#       # Place reshuffle hear if required later. Not need so far     
#        ori_idx = range(len(idx2)+1)
#        new_idx = range(1,pos+1)+[0]+range(pos+1,len(idx2)+1)
#        #print "ori_idx = ",ori_idx
#        #print "new_idx = ",new_idx
#    
#        newmat[:,ori_idx] = newmat[:,new_idx]
#        newmat[ori_idx,:] = newmat[new_idx,:]
#    
        return genCoeffDirect(newmat,CropClass),newmat


def CalculateMatrix(mat,idx1,idx2):
    # Creat reduiced matrix by summing rows and colums according to gathered indeces
    #mtemp=np.sum(mat[:,idx], axis = 1)
    mcol  = np.sum(mat[:,idx1], axis = 1, keepdims = True)
    mkeep = mat[:,idx2]
    mtemp = np.hstack((mcol,mkeep))
    mcol  = np.sum(mtemp[idx1,:], axis = 0, keepdims = True)
    mkeep = mtemp[idx2,:]
    mnew = np.vstack((mcol,mkeep))
    return mnew

def CalculateMatrixMultiple(mat,gatherIdx):
    # Creat reduiced matrix by summing rows and colums according to gathered indeces
    #mtemp=np.sum(mat[:,idx], axis = 1)

    mcol = []
    for idx in gatherIdx:
        mcol.append(np.sum(mat[:,idx], axis = 1, keepdims = True)) 
    mcol = tuple(mcol) 

    mtemp = np.hstack(mcol)
    
    mrow = []
    for idx in gatherIdx:
        mrow.append(np.sum(mtemp[idx,:], axis = 0, keepdims = True))
    mrow = tuple(mrow) 

    mnew = np.vstack(mrow)
    return mnew



def genBinaryStat(cMatrix,CropMask):
	[AllClass,AllClass_prod] = getAllClass(cMatrix)

	mat = confCoordinatesCSV(cMatrix)
	d = defaultdict(list)
 
	for k,v in mat:
		d[k].append(v)

	csv_f = list(d.items())

	mat_sort = gen_confusionMatrix(csv_f,AllClass,AllClass_prod)
        
	nbrGood = mat_sort.trace()
	nbrSample = mat_sort.sum()

        BinaryMatrix = np.zeros((2,2))
        for i,ci in enumerate(CropMask):
            for j,cj in enumerate(CropMask):
                if ci == 0 and cj == 0 : BinaryMatrix[0,0] = BinaryMatrix[0,0] + mat_sort[i][j]
                if ci == 0 and cj == 1 : BinaryMatrix[0,1] = BinaryMatrix[0,1] + mat_sort[i][j]
                if ci == 1 and cj == 0 : BinaryMatrix[1,0] = BinaryMatrix[1,0] + mat_sort[i][j]
                if ci == 1 and cj == 1 : BinaryMatrix[1,1] = BinaryMatrix[1,1] + mat_sort[i][j]

        #print mat_sort
        TN = BinaryMatrix[0,0]
        FN = BinaryMatrix[0,1]
        TP = BinaryMatrix[1,1]
        FP = BinaryMatrix[1,0]
        BinPre = np.array([TN/(TN+FN),TP/(TP+FP)])  
        BinRec = np.array([TN/(TN+FP),TP/(TP+FN)])
	return BinPre,BinRec

def genBinaryStatPerClass(cMatrix,CropMask):
	[AllClass,AllClass_prod] = getAllClass(cMatrix)

        # Order confusion matrix, just in case
	mat = confCoordinatesCSV(cMatrix)
	d = defaultdict(list)
	for k,v in mat:
		d[k].append(v)
	csv_f = list(d.items())
	mat_sort = gen_confusionMatrix(csv_f,AllClass,AllClass_prod)

        a = [] 
        ReVec = []
        PosMRVec = []
        NegMRVec = []
        PrVec = []
        PosFDRVec = []
        NegFDRVec = []
        for i,ci in enumerate(AllClass):
          ri = sum(mat_sort[i,:])    # Total number of reference pixel for the class
          pi = sum(mat_sort[:,i])    # Total number of predicted pixel for the class
          bi = CropMask[i]

          relativeCM = abs(CropMask+bi-1)  # Relaitive crop mask i.e. 1 if class is similar to current class
 
          mat_temp = mat_sort.copy()

          mat_temp[i,i] = 0                         # Avoid double conting of the diagonal term
          PosFN = sum(mat_temp[i,:]*relativeCM)     # Positive Miss pixels
          NegFN = sum(mat_temp[i,:]*(1-relativeCM)) # Negative Miss pixels
          Re = mat_sort[i,i]                        # Recall

          PosFP = sum(mat_temp[:,i]*relativeCM)     # Positive False pixels
          NegFP = sum(mat_temp[:,i]*(1-relativeCM)) # Negative False pixels
          Pr = mat_sort[i,i]                        # Precision

          Re = 100*float(Re)/float(ri)
          PosMR = 100*float(PosFN)/float(ri)
          NegMR = 100*float(NegFN)/float(ri)

          Pr = 100*float(Pr)/float(pi)
          PosFDR = 100*float(PosFP)/float(pi)
          NegFDR = 100*float(NegFP)/float(pi)
  
          ReVec.append(Re)
          PosMRVec.append(PosMR)
          NegMRVec.append(NegMR)
          PrVec.append(Pr)
          PosFDRVec.append(PosFDR)
          NegFDRVec.append(NegFDR)

          #print "ici:",ci, Re, PosMR, NegMR, Re+PosMR+NegMR, Pr, PosFDR, NegFDR, Pr+PosFDR+NegFDR
          #print
          
        return ReVec,PosMRVec,NegMRVec,PrVec,PosFDRVec,NegFDRVec 
  
def genBinaryStatDirect(cMatrix):
	mat = np.loadtxt(cMatrix, delimiter=',', comments='#')

        TN = mat[0,0]
        FN = mat[0,1]
        TP = mat[1,1]
        FP = mat[1,0]
        BinPre = np.array([TN/(TN+FN),TP/(TP+FP)])  
        BinRec = np.array([TN/(TN+FP),TP/(TP+FN)])
	return BinPre,BinRec


def ComputeConfidenceInterval(Data):
	# Degree of freedom
	df = len(Data) - 1

	# Compute standard deviation
	sigma = stdev(Data)
 	if df == 4:
		TDistributionTable = 2.776
	else:
 		if df == 9:
			TDistributionTable = 2.262
		else:
			raise Exception(" Problem")
	return TDistributionTable* (sigma/sqrt(df +1))

class StatisticalStruct(object):
	def __init__(self, Conf, OAmoy, OAint, FSmoy, FSint, PREmoy, PREint, RECmoy, RECint, BinPREmoy, BinPREint, BinRECmoy, BinRECint, dates):
 
                self.Conf = Conf; 
 
		self.OAmoy = OAmoy
		self.OAint = OAint

		self.FSmoy =  FSmoy
		self.FSint =  FSint

                self.PREmoy =  PREmoy
		self.PREint =  PREint

               	self.RECmoy =  RECmoy
		self.RECint =  RECint
	
                self.BinPREmoy =  BinPREmoy
		self.BinPREint =  BinPREint

               	self.BinRECmoy =  BinRECmoy
		self.BinRECint =  BinRECint
		
                self.dates =  dates

class PerClassStatisticalStruct(object):
	def __init__(self, PreMoy,RecMoy,PosMRMoy,NegMRMoy,PosFDRMoy,NegFDRMoy,Dates):
                self.PreMoy = PreMoy
                self.RecMoy = RecMoy
                self.PosMRMoy = PosMRMoy
                self.NegMRMoy = NegMRMoy
                self.PosFDRMoy = PosFDRMoy
                self.NegFDRMoy = NegFDRMoy	
                self.Dates =  Dates

                
def FscoreResults(directory,RFDir,NbDates,Nbtirages,btype,cropmask,impCM,normalized=False):
    Confvec = []
    OAvec = []
    FSvec =[]
    Prevec =[]
    Recvec =[]
    BinPrevec = []
    BinRecvec = []
    for run in range(Nbtirages):
      conf = []
      OA = []
      FS = []
      Pre = []
      Rec = []
      BinPre = [] 
      BinRec = [] 
      t = []
      datesidx= []
      #for j in range(start,NbDates+1,step):
      for j in range(1,NbDates+1):
        iDate = j
        chemin = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_Date%i.csv" %(directory,run,RFDir,btype,iDate)
        #chemin = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_filtered_Date%i.csv" %(directory,run,RFDir,btype,iDate) 
        try:
          #if len(cropmask)>1:
          if btype == "_CM_S2AGRI":
            #BPre, BRec = genBinaryStat(chemin,cropmask)
            BPre, BRec = genBinaryStatDirect(chemin) 
            #cheminCM = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_Date%i.csv" %(directory,run,RFDir,"_CM",iDate)
            #if btype == "_OSORE":
            #    cheminCM = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_Date%i.csv" %(directory,run,RFDir,"_CM_OSORE",iDate)
            #elif btype == "_Fusion_OSORE":
            #    cheminCM = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_Date%i.csv" %(directory,run,RFDir,"_CM_Fusion_OSORE",iDate)
            #BPre, BRec = genBinaryStatDirect(cheminCM) 

            #print btype
            #print chemin 
            chemin = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_Date%i.csv" %(directory,run,RFDir,"_OSORE",iDate)
            #print chemin 
            #print
          else:
            BPre, BRec = genBinaryStat(chemin,cropmask)
            #BPre = np.array([1,1])
            #BRec = np.array([1,1])
          # Exception to deal with imperfect crop mask  
          if impCM:
              #measure = genCoeff(chemin)
              #conf.append(np.loadtxt(chemin,delimiter =',' ,comments='#'))
              measure, mat = genCoeffImperfect(chemin)
              conf.append(mat)
          else:
              measure = genCoeff(chemin,normalized)
              conf.append(np.loadtxt(chemin,delimiter =',' ,comments='#'))
          OA.append(measure.OA)
          FS.append(measure.Fs)
          Pre.append(measure.Pre)
          Rec.append(measure.Rec)
          BinPre.append(BPre)
          BinRec.append(BRec)
          fscore = [ fs[1][:] for fs in FS]
          t.append(iDate)
          datesidx.append(iDate-1)
        except:
          nothing =0

      for i in range(len(np.transpose(FS)[1])):
        labelidx = str(int(np.transpose(FS)[0][i][0]))

      Confvec.append(conf)
      OAvec.append(OA)
      FSvec.append(FS)
      Prevec.append(Pre)
      Recvec.append(Rec)
      BinPrevec.append(BinPre)
      BinRecvec.append(BinRec)
    OAmatrix = np.array(OAvec)
    Conf = np.array(Confvec)

    for i in range(len(OAmatrix.transpose())):
      data = OAmatrix.transpose()[i]
      #data = Conf.transpose()[i]

    AverageConf = np.mean(Conf,axis=0)
    AverageOA = np.mean(OAmatrix,axis=0)
    
    Interval = CI95(OAmatrix)
    FSmatrix = np.array(FSvec)
    Prematrix = np.array(Prevec)
    Recmatrix = np.array(Recvec)
    BinPrematrix = np.array(BinPrevec)
    BinRecmatrix = np.array(BinRecvec)

    FSmoy = np.mean(FSmatrix,axis=0)
    FSInterval = CI95(FSmatrix)
   
    Premoy = np.mean(Prematrix,axis=0)
    PreInterval = CI95(Prematrix)
   
    Recmoy = np.mean(Recmatrix,axis=0)
    RecInterval = CI95(Recmatrix)

    BinPremoy = np.mean(BinPrematrix,axis=0)
    BinPreInterval = CI95(BinPrematrix)
   
    BinRecmoy = np.mean(BinRecmatrix,axis=0)
    BinRecInterval = CI95(BinRecmatrix)
 
    #return [AverageOA,Interval,FSmoy,FSInterval,,,datesidx]
    return StatisticalStruct(AverageConf,AverageOA,Interval,FSmoy,FSInterval,Premoy,PreInterval,Recmoy,RecInterval,BinPremoy,BinPreInterval,BinRecmoy,BinRecInterval,datesidx)

def FscoreResultsGather(directory,RFDir,NbDates,Nbtirages,btype,cropmask,impCM,gatheringSet):
    Confvec = []
    OAvec = []
    FSvec =[]
    Prevec =[]
    Recvec =[]
    BinPrevec = []
    BinRecvec = []
    for run in range(Nbtirages):
      conf = []
      OA = []
      FS = []
      Pre = []
      Rec = []
      BinPre = [] 
      BinRec = [] 
      t = []
      datesidx= []
      #for j in range(start,NbDates+1,step):
      for j in range(1,NbDates+1):
        iDate = j
        chemin = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_Date%i.csv" %(directory,run,RFDir,btype,iDate) 
        #chemin = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_filtered_Date%i.csv" %(directory,run,RFDir,btype,iDate) 
        try:
          #if len(cropmask)>1:
          if btype == "_CM_S2AGRI":
            #BPre, BRec = genBinaryStat(chemin,cropmask)
            BPre, BRec = genBinaryStatDirect(chemin) 
            #cheminCM = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_Date%i.csv" %(directory,run,RFDir,"_CM",iDate)
            #if btype == "_OSORE":
            #    cheminCM = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_Date%i.csv" %(directory,run,RFDir,"_CM_OSORE",iDate)
            #elif btype == "_Fusion_OSORE":
            #    cheminCM = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_Date%i.csv" %(directory,run,RFDir,"_CM_Fusion_OSORE",iDate)
            #BPre, BRec = genBinaryStatDirect(cheminCM) 

            #print btype
            #print chemin 
            chemin = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_Date%i.csv" %(directory,run,RFDir,"_OSORE",iDate)
            #print chemin 
            #print
          else:
            BPre, BRec = genBinaryStat(chemin,cropmask)
            #BPre = np.array([1,1])
            #BRec = np.array([1,1])
          # Exception to deal with imperfect crop mask  
          if impCM:
              #measure = genCoeff(chemin)
              #conf.append(np.loadtxt(chemin,delimiter =',' ,comments='#'))
              measure, mat = genCoeffImperfect(chemin)
              conf.append(mat)
          else:
              measure, mat = genCoeffGather(chemin,gatheringSet)
              conf.append(mat)
 
          OA.append(measure.OA)
          FS.append(measure.Fs)
          Pre.append(measure.Pre)
          Rec.append(measure.Rec)
          BinPre.append(BPre)
          BinRec.append(BRec)
          fscore = [ fs[1][:] for fs in FS]
          t.append(iDate)
          datesidx.append(iDate-1)
        except:
          nothing =0
      for i in range(len(np.transpose(FS)[1])):
        labelidx = str(int(np.transpose(FS)[0][i][0]))

      Confvec.append(conf)
      OAvec.append(OA)
      FSvec.append(FS)
      Prevec.append(Pre)
      Recvec.append(Rec)
      BinPrevec.append(BinPre)
      BinRecvec.append(BinRec)
     
    OAmatrix = np.array(OAvec)
    Conf = np.array(Confvec)

    for i in range(len(OAmatrix.transpose())):
      data = OAmatrix.transpose()[i]
      #data = Conf.transpose()[i]

    AverageConf = np.mean(Conf,axis=0)
    AverageOA = np.mean(OAmatrix,axis=0)
    
    Interval = CI95(OAmatrix)
    FSmatrix = np.array(FSvec)
    Prematrix = np.array(Prevec)
    Recmatrix = np.array(Recvec)
    BinPrematrix = np.array(BinPrevec)
    BinRecmatrix = np.array(BinRecvec)

    FSmoy = np.mean(FSmatrix,axis=0)
    FSInterval = CI95(FSmatrix)
   
    Premoy = np.mean(Prematrix,axis=0)
    PreInterval = CI95(Prematrix)
   
    Recmoy = np.mean(Recmatrix,axis=0)
    RecInterval = CI95(Recmatrix)

    BinPremoy = np.mean(BinPrematrix,axis=0)
    BinPreInterval = CI95(BinPrematrix)
   
    BinRecmoy = np.mean(BinRecmatrix,axis=0)
    BinRecInterval = CI95(BinRecmatrix)
 
    #return [AverageOA,Interval,FSmoy,FSInterval,,,datesidx]
    return StatisticalStruct(AverageConf,AverageOA,Interval,FSmoy,FSInterval,Premoy,PreInterval,Recmoy,RecInterval,BinPremoy,BinPreInterval,BinRecmoy,BinRecInterval,datesidx)


def PerClassResults(directory,RFDir,NbDates,Nbtirages,btype,cropmask):
    PreVec =[]
    RecVec =[]
    PosMissRateVec = []
    NegMissRateVec = []
    PosFalseDRateVec = []
    NegFalseDRateVec = []
    for run in range(Nbtirages):
      Rec = []
      Pre = []
      PosMissRate = []
      NegMissRate = []
      PosFalseDRate = []
      NegFalseDRate = []

      t = []
      datesidx= []
      #for j in range(start,NbDates+1,step):
      for j in range(1,NbDates+1):
      #for j in [2,6,10]: # For debug
        iDate = j
        chemin = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_Date%i.csv" %(directory,run,RFDir,btype,iDate) 
        #chemin = "%sRun_%i.dir/%s/ConfusionsMatrices/ConfusionsMatrix%s_filtered_Date%i.csv" %(directory,run,RFDir,btype,iDate) 
        try:
          #BPre, BRec = genBinaryStat(chemin,cropmask)
          Re,PosMR,NegMR,Pr,PosFDR,NegFDR = genBinaryStatPerClass(chemin,cropmask)
          #conf.append(np.loadtxt(chemin,delimiter =',' ,comments='#'))
          
          Rec.append(Re)
          Pre.append(Pr)
          PosMissRate.append(PosMR)
          NegMissRate.append(NegMR)
          PosFalseDRate.append(PosFDR)
          NegFalseDRate.append(NegFDR) 
          
          #fscore = [ fs[1][:] for fs in FS]
          t.append(iDate)
          datesidx.append(iDate-1)
        except:
          nothing = 0 
      #for i in range(len(np.transpose(FS)[1])):
      #  labelidx = str(int(np.transpose(FS)[0][i][0]))

      PreVec.append(Rec)
      RecVec.append(Pre)
      PosMissRateVec.append(PosMissRate)
      NegMissRateVec.append(NegMissRate)
      PosFalseDRateVec.append(PosFalseDRate)
      NegFalseDRateVec.append(NegFalseDRate)


    
   #OAmatrix = np.array(OAvec)
   #Conf = np.array(Confvec)

   # for i in range(len(OAmatrix.transpose())):
   #   data = OAmatrix.transpose()[i]
   #   #data = Conf.transpose()[i]

    PreMatrix = np.array(PreVec)
    RecMatrix = np.array(RecVec)    
    PosMRMatrix = np.array(PosMissRateVec)
    NegMRMatrix = np.array(NegMissRateVec)
    PosFDRMatrix = np.array(PosFalseDRateVec)
    NegFDRMatrix = np.array(NegFalseDRateVec)
    PreMoy = np.mean(PreMatrix,axis=0)
    RecMoy = np.mean(RecMatrix,axis=0)
    PosMRMoy = np.mean(PosMRMatrix,axis=0)
    NegMRMoy = np.mean(NegMRMatrix,axis=0)
    PosFDRMoy = np.mean(PosFDRMatrix,axis=0)
    NegFDRMoy = np.mean(NegFDRMatrix,axis=0)

 
    #FSmoy = np.mean(FSmatrix,axis=0)
    #FSInterval = CI95(FSmatrix)
   
    #Premoy = np.mean(Prematrix,axis=0)
    #PreInterval = CI95(Prematrix)
   
    #Recmoy = np.mean(Recmatrix,axis=0)
    #RecInterval = CI95(Recmatrix)

    #BinPremoy = np.mean(BinPrematrix,axis=0)
    #BinPreInterval = CI95(BinPrematrix)
   
    #BinRecmoy = np.mean(BinRecmatrix,axis=0)
    #BinRecInterval = CI95(BinRecmatrix)
 
    return PerClassStatisticalStruct(PreMoy,RecMoy,PosMRMoy,NegMRMoy,PosFDRMoy,NegFDRMoy,datesidx)

def PerClassTableau(cname,data1,data2,data3,doy,cap,mark):

    tf = "|l||c|c|c||c|c|c||c|c|c|"
    mc = ["",doy[0],doy[1],doy[2]]

    nd = data1.shape[0]
    nc = data1.shape[1]
    print nd,nc
    print data2.shape
    print data3.shape


    # x = [fl[0][d][i][1],fl[1][d][i][1],fl[2][d][i][1]]
    # idx = x.index(max(x))
    # if idx == 0: fl[0][d][i][1] = "\\bf{%.2f}"%(fl[0][d][i][1]) 
    # elif idx == 1: fl[1][d][i][1] = "\\bf{%.2f}"%(fl[1][d][i][1]) 
    # elif idx == 2: fl[2][d][i][1] = "\\bf{%.2f}"%(fl[2][d][i][1]) 



    alldata = [["Class","Opt","Rad","Fus","Opt","Rad","Fus","Opt","Rad","Fus"]]
    for c in range(nc):
        line = [cname[c],data1[0][c],data2[0][c],data3[0][c],data1[1][c],data2[1][c],data3[1][c],data1[2][c],data2[2][c],data3[2][c]]
        for i in [1,4,7]:
            if  (mark=="Max"):idx = line.index(max(line[i:i+3]))
            elif(mark=="Min"):idx = line.index(min(line[i:i+3]))
            else:
                print("mark should be max or min")
                quit()
            line[idx] = "\\bf{%.2f}"%(line[idx])
        if(line[0] != "Forest"): alldata.append(line)

    #for i in l
    doc2.tableau([alldata],tf, multicol = mc, size = [1,3,3,3],caption = cap + " (%simum are indicated with bold font)."%(mark))
    



def printFscore(fscoreRADtmp,RadOA,RadInt,fscoreOPTtmp,OptOA,OptInt,fscoreFUStmp,FusOA,FusInt,RadDates,OptDates,DoY,classname,document,captionstring):
        
    tf = "|l||r|r|r||r|r|r||r|r|r|"

    fscoreRAD = np.zeros((fscoreRADtmp.shape[0],fscoreRADtmp.shape[1]+2,2))
    fscoreOPT = np.zeros((fscoreOPTtmp.shape[0],fscoreOPTtmp.shape[1]+2,2))
    fscoreFUS = np.zeros((fscoreFUStmp.shape[0],fscoreFUStmp.shape[1]+2,2))

    for t in range(len(fscoreRADtmp)):
      fscoreRAD[t] = np.append(fscoreRADtmp[t], [[-1,RadOA[t]],[-2,RadInt[t]]],axis = 0)
      fscoreOPT[t] = np.append(fscoreOPTtmp[t], [[-1,OptOA[t]],[-2,OptInt[t]]],axis = 0)
      fscoreFUS[t] = np.append(fscoreFUStmp[t], [[-1,FusOA[t]],[-2,FusInt[t]]],axis = 0)

    fscoreRAD = 100.0*fscoreRAD
    fscoreOPT = 100.0*fscoreOPT
    fscoreFUS = 100.0*fscoreFUS

    fl=[fscoreRAD[RadDates].tolist(),fscoreOPT[OptDates].tolist(),fscoreFUS[RadDates].tolist()]
    print len(fl[0][0])
    finaldata = []
    i=0

    newblock = []
    newline = ["Class",0,0,0,0,0,0,0,0,0]
    newblock.append(newline)
    finaldata.append(newblock) 
    for f in (fl[0][0]):
        newblock = []
        
        for d in [0,1,2]:
           x = [fl[0][d][i][1],fl[1][d][i][1],fl[2][d][i][1]]
           idx = x.index(max(x))
           if idx == 0: fl[0][d][i][1] = "\\bf{%.2f}"%(fl[0][d][i][1]) 
           elif idx == 1: fl[1][d][i][1] = "\\bf{%.2f}"%(fl[1][d][i][1]) 
           elif idx == 2: fl[2][d][i][1] = "\\bf{%.2f}"%(fl[2][d][i][1]) 


        newline = [classname[int(f[0])/100],fl[0][0][i][1],fl[1][0][i][1],fl[2][0][i][1],fl[0][1][i][1],fl[1][1][i][1],fl[2][1][i][1],fl[0][2][i][1],fl[1][2][i][1],fl[2][2][i][1]]
        newblock.append(newline)
        finaldata.append(newblock)
        i = i + 1
 
    multicol = ["", DoY[0], DoY[1], DoY[2]]
    finaldata[0][0][1] = "Rad"; finaldata[0][0][4] = "Rad";finaldata[0][0][7] = "Rad"
    finaldata[0][0][2] = "Opt"; finaldata[0][0][5] = "Opt"; finaldata[0][0][8] = "Opt"
    finaldata[0][0][3] = "Fus"; finaldata[0][0][6] = "Fus"; finaldata[0][0][9] = "Fus"

    document.tableau(finaldata,tf, multicol = multicol, size = [1,3,3,3],caption = captionstring)

def printFscoreNew(RadStat,OptStat,FusStat,RadDates,OptDates,DoY,classname,document,captionstring):
   
    fscoreRADtmp  = RadStat.FSmoy
    RadOA         = RadStat.OAmoy
    RadInt        = RadStat.OAint

    fscoreOPTtmp  = OptStat.FSmoy
    OptOA         = OptStat.OAmoy
    OptInt        = OptStat.OAint

    fscoreFUStmp  = FusStat.FSmoy
    FusOA         = FusStat.OAmoy
    FusInt        = FusStat.OAint

    tf = "|l||r|r|r||r|r|r||r|r|r|"

    fscoreRAD = np.zeros((fscoreRADtmp.shape[0],fscoreRADtmp.shape[1]+2,2))
    fscoreOPT = np.zeros((fscoreOPTtmp.shape[0],fscoreOPTtmp.shape[1]+2,2))
    fscoreFUS = np.zeros((fscoreFUStmp.shape[0],fscoreFUStmp.shape[1]+2,2))

    for t in range(len(fscoreRADtmp)):
      fscoreRAD[t] = np.append(fscoreRADtmp[t], [[-1,RadOA[t]],[-2,RadInt[t]]],axis = 0)
      fscoreOPT[t] = np.append(fscoreOPTtmp[t], [[-1,OptOA[t]],[-2,OptInt[t]]],axis = 0)
      fscoreFUS[t] = np.append(fscoreFUStmp[t], [[-1,FusOA[t]],[-2,FusInt[t]]],axis = 0)

    fscoreRAD = 100.0*fscoreRAD
    fscoreOPT = 100.0*fscoreOPT
    fscoreFUS = 100.0*fscoreFUS

    fl=[fscoreRAD[RadDates].tolist(),fscoreOPT[OptDates].tolist(),fscoreFUS[RadDates].tolist()]
    print len(fl[0][0])
    finaldata = []
    i=0

    newblock = []
    newline = ["Class",0,0,0,0,0,0,0,0,0]
    newblock.append(newline)
    finaldata.append(newblock) 
    for f in (fl[0][0]):
        newblock = []
        
        for d in [0,1,2]:
           x = [fl[0][d][i][1],fl[1][d][i][1],fl[2][d][i][1]]
           idx = x.index(max(x))
           if idx == 0: fl[0][d][i][1] = "\\bf{%.2f}"%(fl[0][d][i][1]) 
           elif idx == 1: fl[1][d][i][1] = "\\bf{%.2f}"%(fl[1][d][i][1]) 
           elif idx == 2: fl[2][d][i][1] = "\\bf{%.2f}"%(fl[2][d][i][1]) 


        newline = [classname[int(f[0])/100],fl[0][0][i][1],fl[1][0][i][1],fl[2][0][i][1],fl[0][1][i][1],fl[1][1][i][1],fl[2][1][i][1],fl[0][2][i][1],fl[1][2][i][1],fl[2][2][i][1]]
        newblock.append(newline)
        finaldata.append(newblock)
        i = i + 1
 
    multicol = ["", DoY[0], DoY[1], DoY[2]]
    finaldata[0][0][1] = "Rad"; finaldata[0][0][4] = "Rad";finaldata[0][0][7] = "Rad"
    finaldata[0][0][2] = "Opt"; finaldata[0][0][5] = "Opt"; finaldata[0][0][8] = "Opt"
    finaldata[0][0][3] = "Fus"; finaldata[0][0][6] = "Fus"; finaldata[0][0][9] = "Fus"

 

    document.tableau(finaldata,tf, multicol = multicol, size = [1,3,3,3],caption = captionstring)



def classcount(x):
   classes = []
   nbr = [] 
   for i in range(len(x)):
      cl = int(x[i])
      if cl not in classes:
         classes.append(cl)  
         nbr.append(1)
      else:
         nbr[classes.index(cl)] += 1
   z = np.array((classes,nbr))
   z = z.transpose()
   z = np.array(sorted(z, key=lambda row: row[0]))
   return z

def printline(l):
   for i in range(len(l.transpose()[0])):
      print l.transpose()[0][i],"\t",
   print
   for i in range(len(l.transpose()[0])):
      print l.transpose()[1][i],"\t",
   print

def GetParcellCount(shapeFile):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapeFile, 0)
    layer = dataSource.GetLayer()
    
    classes = []
    ParcelleNbr = []
    for feature in layer:
        idCode = int(feature.GetField("CODE"))
        if idCode not in classes:
            classes.append(idCode)
            ParcelleNbr.append(1)
        else:
            ParcelleNbr[classes.index(idCode)] += 1
    z = np.array((classes,ParcelleNbr))
    z = z.transpose()
    z = np.array(sorted(z, key=lambda row: row[0]))
    classes = z.transpose()[0].tolist()
    ParcelleNbr = z.transpose()[1].tolist()
    return classes,ParcelleNbr


def PrintConfMat(RFDir,btype,captionstring,document):
    ConfMat = 0
    for i in range(Nbtirages):
        chemin = "%sRun_%i/%s/ConfuMatri/%sConfuMatrix_Date%i.csv" %(directory,i,RFDir,btype,33)
        #chemin = "%sRun_%i/%s/ConfuMatri/%sConfuMatrix_filtered_Date%i.csv" %(directory,i,RFDir,btype,33)
        ConfMat = np.add(ConfMat, loadtxt(chemin, delimiter = ","))
        AllClass,AllClass_prod = getAllClass(chemin)
    
    ClassName = [ classname[key] for key in AllClass]
   
    # To sort the matrix a the nice way Crop/Nocrop
    if len(ClassName) > 8:
        tmp = ClassName[0]
        ClassName[0] = ClassName[10]
        ClassName[10] = tmp
        tmp = AllClass[0]
        AllClass[0] = AllClass[10]
        AllClass[10] = tmp
        

    ClassName = [""] +  ClassName
    mat = (ConfMat/Nbtirages).astype(int).tolist()

    if len(ClassName) > 8:
        newmat = mat
        newmat[0] = mat[10]
        newmat[10] = mat[0]
    else:
       newmat = mat
        

    
    for i in range(len(newmat)):
        newmat[i] = [classname[AllClass[i]]] + newmat[i]
    
    document.matrix([[ClassName],newmat],"|l|r|r|r|r|r|r|r|r|r|r|r|r|r|r|r|r|",caption = captionstring)


def printmatrix(ClassNames,mat,cap,doc):
    ClassNames = ["Class"] + ClassNames
    listmat = mat.astype(int).tolist()
    newmat = listmat
    tf = "|l|"
    for i in range(len(newmat)):
      tf = tf + "r|"
      newmat[i] = [ClassNames[i+1]] + newmat [i] 
    
    finalmat = [[ClassNames],newmat]
    
    doc.matrix(finalmat,tf,caption=cap)

    
if __name__ == '__main__':
    
    ############################################################################
    
    # General Parameters #
    
    # Band labels
    bandstype = {
    "":"Radar",
    "_filtered":"Radar",
    "_ROM":"Radar Mix. Orb",
    "_ROM_filtered":"Radar Mix. Orb",
    "_WSC":"Radar",
    "_IMP":"Radar",
    "_WSC_IMP":"Radar",
    "_CROP":"Optical",
    "_CROP_WSC":"Optical",
    "_CROP_IMP":"Optical",
    "_CROP_WSC_IMP":"Optical",
    "_OSORE":"Optical+OSO+RedEdge",
    "_OSORE_filtered":"Optical+OSO+RedEdge",
    "_OSORE_IMP":"Optical+OSO+RedEdge",
    "_OSORE_WSC":"Optical+OSO+RedEdge",
    "_OSORE_WSC_IMP":"Optical+OSO+RedEdge",
    "_CROP_OSORE":"Optical+OSO+RedEdge",
    "_CROP_OSORE_WSC":"Optical+OSO+RedEdge",
    "_CM_S2AGRI":"S2AGRI",
    "_OSORE_1e4":"OSORE 1e4",
    "_PCA1":"PCA 1",
    "_PCAk4":"PC4 k=4",
    "_PCAk5":"PCA k=5",
    "_PCA61":"PCA 61",
    "_TEST3":"OSORE TEST3",
    "_Fusion":"Fusion",
    "_Fusion_WSC":"Fusion",
    "_Fusion_WSC_IMP":"Fusion",
    "_Fusion_OSORE":"Fusion",
    "_Fusion_OSORE_filtered":"Fusion",
    "_Fusion_OSORE_IMP":"Fusion",
    "_Fusion_OSORE_WSC":"Fusion",
    "_Fusion_OSORE_WSC_IMP":"Fusion",
    "_CROP_Fusion":"Fusion",
    "_CROP_Fusion_WSC":"Fusion",
    "_CROP_Fusion_OSORE":"Fusion",
    "_CROP_Fusion_OSORE_WSC":"Fusion",
    "_S2AGRI":"S2Agri"
    }
    
    # Input parameter
    if len(sys.argv) != 8:
        print "USAGE"
        print "./PlotFscoreAll.py workdir NbRun type title RadDatesFile OptDatesFile ispngout"
        print "type 0  = All Classes"
        print "     1  = Prefect Crop Type"
        print "     2  = Imperfect Crop Type"
        quit()
    
    workdir = sys.argv[1]
    Nbtirages = int(sys.argv[2])
    croptype = int(sys.argv[3])
    title = sys.argv[4]
    RadDatesFile = sys.argv[5]
    OptDatesFile = sys.argv[6]
    png = int(sys.argv[7])
    
    
    classfile = workdir + "/WorkFiles/Classes.csv" 
    gatherfile = workdir + "/WorkFiles/Classes_Gathered.csv" 
    binaryfile = workdir + "/WorkFiles/Binary.txt" 
    colorfile = workdir + "/WorkFiles/Legend.csv" 
    outputdir = workdir + "/Statistics/" 
    outputdir = "" 
    
    # Construct Class name dictionary 
    with open(classfile, "r") as ins:
        classname = {}
        for line in ins:
            couple = line.split(":")
            classname[int(couple[1])] = couple[0]
        classname[0]="NoCrop"
    classname[-1] = "OA"
    classname[-2] = "95 CI"
    
    # Construct Class gathering dictionary 
    with open(gatherfile, "r") as ins:
        gatheringSet = {}
        for line in ins:
            couple = line.split(":")
            gatheringSet[int(couple[0])] = int(couple[1])
        gatheringSet[0]=0
    
    # Construct Color dictionary
    with open(colorfile, "r") as ins:
        colordic = {0:"#000000"}
        for line in ins:
            couple = line.split(":")
            colordic[int(couple[0])] = couple[1][:-1] # [:-1] to remove the "\n" at the end of the line 
    
    # Fix cropmask "mask"
    if croptype==0:
      cropmask = np.loadtxt(binaryfile)
      # filtered
      cropmask = cropmask[:-1].copy()
    
    else: cropmask = np.array([1])
    
    # Path and directory settings
    RFDirOpt = "Classifications"
    RFDirRad = "Classifications"
    optical_dates_file = workdir + "/WorkFiles/" + OptDatesFile
    radar_dates_file = workdir + "/WorkFiles/" + RadDatesFile
    FusDirectory = workdir + "/AllTiles/Fusion/"
    S2AGRIDirectory = workdir + "/AllTiles/RadarOptical/"
    
    # Optical Parameters #
    OpticalTest = True
    OptDirectory = workdir + "/AllTiles/Optical/"
    OptDate,OptDoYlabel = GetDates(optical_dates_file)
    NbOptDates = len(OptDate)
    OptDate = np.array(OptDate)
    DoYDic = {}
    for i,d in enumerate(OptDate):
        DoYDic[d]=OptDoYlabel[i]
    
    # Radar Parameters #
    RadarTest = True
    RadDirectory = workdir + "/AllTiles/Radar/"
    RadDate,RadDoYlabel = GetDates(radar_dates_file)
    NbRadDates = len(RadDate)
    RadDate = np.array(RadDate)
    for i,d in enumerate(RadDate):
        DoYDic[d]=RadDoYlabel[i]
    
    # Fusion Parameters #
    FusionTest = True
    
    # Set variables for "Crop type after crop mask"
    TypeName=["","-CROP","CTAfterCM"]
    if (png==1):os.system("mkdir -p fig")
    
    imperfectCM = False
    # Set input and output files variables
    if(croptype==0):
        if(title=="Spain2017"):
            # Spain special case
            RadBtypeList = ["_WSC"]
            OptBtypeList = ["_OSORE_WSC"]
            FusBtypeList = ["_Fusion_OSORE_WSC"]
        else:
            RadBtypeList = ["_filtered"]
            OptBtypeList = ["_OSORE_filtered"]
            FusBtypeList = ["_Fusion_OSORE_filtered"]
        PlotTitle =  title + ": All Classes Fscore: "
        FscorePDF = outputdir + "Plot-Fscore-AllClasses.pdf"
        ConfusionsPDF = outputdir + "Plot-Confusions-AllClasses.pdf"
        TablePDF = outputdir + "Table-Fscore-AllClasses"
    
    elif(croptype==1):
        if(title=="Spain2017"):
            # Spain special case
            RadBtypeList = ["_CROP_WSC"]
            OptBtypeList = ["_CROP_OSORE_WSC"]
            FusBtypeList = ["_CROP_Fusion_OSORE_WSC"]
        else:
            s2agriBtypeList = ["_S2AGRI"]
            RadBtypeList = ["_CROP"]
            OptBtypeList = ["_CROP_OSORE"] 
            FusBtypeList = ["_CROP_Fusion_OSORE"]
        PlotTitle =  title + ": Crop Type Fscore: "
        FscorePDF = outputdir + "Plot-Fscore-CropType.pdf"
        ConfusionsPDF = outputdir + "Plot-Confusions-CropType.pdf"
        TablePDF = outputdir + "Table-Fscore-CropType"
        TablePDFgather = outputdir + "Table-Fscore-CropType-Gather"
    
    elif(croptype==2):
        imperfectCM = True
        if(title=="Spain2017"):
            # Replacement of Winter and Spring cererals for Spain 2017 #
            #RadBtypeList = ["_WSC_IMP"]
            #OptBtypeList = ["_OSORE_WSC_IMP"]
            #FusBtypeList = ["_Fusion_OSORE_WSC_IMP"]
            RadBtypeList = ["_WSC"]
            OptBtypeList = ["_OSORE_WSC"]
            FusBtypeList = ["_Fusion_OSORE_WSC"]
        else:
            #RadBtypeList = ["_IMP2"]
            #OptBtypeList = ["_OSORE_IMP2"]
            #FusBtypeList = ["_Fusion_OSORE_IMP2"]
            RadBtypeList = [""]
            OptBtypeList = ["_OSORE"]
            FusBtypeList = ["_Fusion_OSORE"]
        PlotTitle =  title + ": Crop type after crop mask Fscore: "
        FscorePDF = outputdir + "Plot-Fscore-CTafterCM2.pdf"
        ConfusionsPDF = outputdir + "Plot-Confusions-CTafterCM2.pdf"
        TablePDF = outputdir + "Table-Fscore-CTafterCM2"
    
    # Loading of optical data in order to get the number of class #
    Stat = FscoreResults(OptDirectory,RFDirOpt,NbOptDates,Nbtirages,OptBtypeList[0],cropmask,imperfectCM)
    OptFscore = Stat.FSmoy
    NbClasses= len(OptFscore[0])
    
    # Open pdf files and set color range#
    if(png == 0): pdf = matplotlib.backends.backend_pdf.PdfPages(FscorePDF)
    #col = ["#AAAAAA","#0000FF","#FFAA00","#005555","#CC0000","#00CC00",]
    col = ["#0000FF","#FFAA00","#CC0000","#00CC00",]
    # Blue for optical
    # Red or Radar 
    # Green for Fusion
    xl = 0.9
    yl = 0.3
    
    # Set measure dictionaries
    OptStatList = {}
    OptPerClassStatList = {}
    RadStatList = {}
    RadPerClassStatList = {}
    FusStatList = {}
    FusPerClassStatList = {}
    for btype in OptBtypeList:
        print "Get Optical Measures"
        OptStatList[btype] = FscoreResults(OptDirectory,RFDirOpt,NbOptDates,Nbtirages,btype,cropmask,imperfectCM)
        #if(croptype == 0): OptPerClassStatList[btype] = PerClassResults(OptDirectory,RFDirOpt,NbOptDates,Nbtirages,btype,cropmask)
    
    for btype in RadBtypeList:
        print "Get Radar Measures"
        print RadDirectory,RFDirRad,NbRadDates,Nbtirages,btype,cropmask,imperfectCM
        RadStatList[btype] = FscoreResults(RadDirectory,RFDirRad,NbRadDates,Nbtirages,btype,cropmask,imperfectCM)
        #if(croptype == 0): RadPerClassStatList[btype] = PerClassResults(RadDirectory,RFDirRad,NbRadDates,Nbtirages,btype,cropmask)
    
    for btype in FusBtypeList:
        print "Get Fusion Measures"
        FusStatList[btype] = FscoreResults(FusDirectory,RFDirOpt,NbOptDates,Nbtirages,btype,cropmask,imperfectCM)
        #if(croptype == 0): FusPerClassStatList[btype] = PerClassResults(FusDirectory,RFDirOpt,NbOptDates,Nbtirages,btype,cropmask)
    
    # Create tables per class Crop Mask statistics
    #if(croptype == 0):
    
    #    Dates = OptPerClassStatList["_OSORE"].Dates
    #    #Dates = OptPerClassStatList["_OSORE_WSC"].Dates
    #    print OptDoYlabel
    #    print  Dates
    #    Idx = [2,4,6]
    #    DoY = [OptDoYlabel[Dates[Idx[0]]],OptDoYlabel[Dates[Idx[1]]],OptDoYlabel[Dates[Idx[2]]]]
    #    print DoY
    #    print FusStatList["_Fusion_OSORE"].BinRECmoy[Idx,0],FusStatList["_Fusion_OSORE"].BinPREmoy[Idx,0]
    #    print FusStatList["_Fusion_OSORE"].BinRECmoy[Idx,1],FusStatList["_Fusion_OSORE"].BinPREmoy[Idx,1]
    #    #print FusStatList["_Fusion_OSORE_WSC"].BinRECmoy[Idx,0],FusStatList["_Fusion_OSORE_WSC"].BinPREmoy[Idx,0]
    #    #print FusStatList["_Fusion_OSORE_WSC"].BinRECmoy[Idx,1],FusStatList["_Fusion_OSORE_WSC"].BinPREmoy[Idx,1]
    #
    #    quit()
    #
    #    #PreMoy,RecMoy,PosMRMoy,NegMRMoy,PosFDRMoy,NegFDRMoy,datesidx
    #    doc2 = latex.latexdoc("PerClassStatistics")
    #    doc2.title("Per Class Crop Mask Statistics")
    #
    ##    Idx = [2,4,6]
    #    #RadIdx = [2,4,6]
    #    Dates = OptPerClassStatList["_OSORE"].Dates
    #    DoY = [OptDoYlabel[Dates[Idx[0]]],OptDoYlabel[Dates[Idx[1]]],OptDoYlabel[Dates[Idx[2]]]]
    #    cname = [classname[int(OptFscore[0][c][0])] for c in range(len(OptFscore[0]))]
    #    doc2.write("\\section{Formal definition and properties}")
    #    doc2.write("\\subsection{Definitions of statistical quantities}")
    #    doc2.write("""\\begin{itemize}
    #                   \\item[$\\bullet$] $P_{c_i}$ is the precision per class defined as
    #                   \\begin{equation}
    #                   P_{c_i}=  \\frac{ \\displaystyle  \\sum_{ \\forall x ~ | ~ \\substack{ c_p =  c_i   \\\\ c_r = c_i }} x } { \\displaystyle  \\sum_{ \\forall x ~ \\in  ~ c_p } x}
    #                   \\end{equation}.
    #                   \\item[$\\bullet$] $R_{c_i}$ is the Recall per class defined as
    #
    #                   \\begin{equation}
    #                   R_{c_i}=  \\frac{ \\displaystyle  \\sum_{ \\forall x ~ | ~ \\substack{ c_p =  c_i   \\\\ c_r = c_i }} x } { \\displaystyle  \\sum_{ \\forall x ~ \\in  ~ c_r } x}
    #                  \\end{equation}.
    #                  \\item[$\\bullet$] $F_{c_i}^{+}$ is the {\\em Positive False Discovery Rate} defined as
    #                  \\begin{equation}
    #    F_{c_i}^{-}=  \\frac{ \\displaystyle  \\sum_{ \\forall x ~ | ~ \\substack{ c_p =  c_j   \\\\ c_r = c_i  \\\\ c_j \\in -  }} x } { \\displaystyle  \\sum_{ \\forall x ~ \\in  ~ c_p } x}
    #    \end{equation}
    #                  \\item[$\\bullet$] $F_{c_i}^{-}$ is the {\\em Negative False Discovery Rate} defined as
    #                  \\begin{equation}
    #    F_{c_i}^{+}=  \\frac{ \\displaystyle  \\sum_{ \\forall x ~ | ~ \\substack{ c_p =  c_j   \\\\ c_r = c_i  \\\\ c_j \\in +  }} x } { \\displaystyle  \\sum_{ \\forall x ~ \\in  ~ c_p } x}
    #    \end{equation}
    #                  \\item[$\\bullet$] $M_{c_i}^{+}$ is the {\\em Positive Miss Rate} defined as
    #                  \\begin{equation}
    #                  M_{c_i}^{+}=  \\frac{ \\displaystyle  \\sum_{ \\forall x ~ | ~ \\substack{ c_p /=  c_r   \\\\ c_r = c_i  \\\\ c_r \\in +  }} x } { \\displaystyle  \\sum_{ \\forall x \in c_r}   x}
    #                  \end{equation}
    #                  \\item[$\\bullet$] $M_{c_i}^{-}$ is the {\\em Negative Miss Rate} defined as
    #                  \\begin{equation}
    #                  M_{c_i}^{-}=  \\frac{ \\displaystyle  \\sum_{ \\forall x ~ | ~ \\substack{ c_p /=  c_r   \\\\ c_r = c_i  \\\\ c_r \\in -  }} x } { \\displaystyle  \\sum_{ \\forall x \in c_r}   x}
    #    \end{equation}
    #                  \\end{itemize}
    #                 """)
    #    doc2.write("\\subsection{Properties}")
    #    doc2.write("""\\begin{itemize}
    #                 \\item $P_{c_i} +F_{c_i}^{+} + F_{c_i}^{-} = 100$
    #                 \\item $R_{c_i} +M_{c_i}^{+} + M_{c_i}^{-} = 100$
    #                 \\end{itemize}
    #    """)
    #     
    #    doc2.write("\\newpage")
    #    doc2.write("\\section{Intuitive definition thanks to the confusion matrix}")
    #    doc2.write("""Let's consider a confusion matrix $C$, which is the matrix where the element $C_{ij}$ is the number of reference pixels of the class $i$ that is predicted to belongs to the class $j$. Summing the elements on a given row $i$ gives the total number of reference pixels of the class $i$:
    #    $$ \\sum_{j}C_{ij} = r_j = \\textrm{Number of reference pixels of the class } i.$$
    #    Similarly, summing the elements on a given column $j$ gives the total number of pixels predicted to belong to the class $j$:
    #    $$ \\sum_{i}C_{ij} = p_j = \\textrm{Number of predicted pixels of the class } j.$$
    #    Each row and column can be decomposed in a way that enphasis the different kind of confution between pixel. By introducing the following counting of the pixels,
    #    \\begin{itemize}
    #    \\item The diagonal term $C_{ii}$ as the number of correclty classi.
    #    \\item The number of in correctly classified pixels.
    #    \\end{itemize}
    #
    # Following this decomposition, we have:
    #    $$ r_i = C_{ii} + FN^{+}_i + FN^{-}_i$$
    #    $$ p_j = C_{jj} + FP^{+}_j + FP^{-}_j$$
    # 
    #By renormalising each term, we got the following expression for the statistical quantities of interrest:
    #    $$ R_i = C_{ii}/r_i$$
    #    $$ M^{+}_i = FN^{+}_i/r_i$$
    #    $$ M^{-}_i = FN^{-}_i/r_i$$
    #    $$ P_j = C_{jj}/p_j$$
    #    $$ F^{+}_j = FP^{+}_j/p_j$$
    #    $$ F^{-}_j = FP^{-}_j/p_j$$
    #
    #""")
    #    
    # 
    #    # Recal per class
    #    doc2.write("\\newpage")
    #    doc2.write("\\section{Recall per class}")
    #    doc2.write("\\subsection{Definition}")
    #    doc2.write("Recall per class: Proportion of pixels correctly classify among the reference data.")
    #    doc2.write("\\subsection{Results}")
    #    data1 = OptPerClassStatList["_OSORE"].RecMoy[Idx]
    #    data2 = RadPerClassStatList[""].RecMoy[Idx]
    #    data3 = FusPerClassStatList["_Fusion_OSORE"].RecMoy[Idx]
    #    PerClassTableau(cname,data1,data2,data3,DoY, title +": Recal per class averaged over 10 runs","Max")
    # 
    #    # Positive miss rate per class
    #    doc2.write("\\newpage")
    #    doc2.write("\\section{Positive miss rate}")
    #    doc2.write("\\subsection{Definition}")
    #    doc2.write("Positive miss rate: Proportion of pixels lost to an other class of the same crop/nocrop kind.")
    #    doc2.write("\\subsection{Results}")
    #    data1 = OptPerClassStatList["_OSORE"].PosMRMoy[Idx]
    #    data2 = RadPerClassStatList[""].PosMRMoy[Idx]
    #    data3 = FusPerClassStatList["_Fusion_OSORE"].PosMRMoy[Idx]
    #    PerClassTableau(cname,data1,data2,data3,DoY, title +": Positive miss rate averaged over 10 runs","Min")
    #
    #    # Negative miss rate per class
    #    doc2.write("\\newpage")
    #    doc2.write("\\section{Negative miss rate}")
    #    doc2.write("\\subsection{Definition}")
    #    doc2.write("Negative miss rate: Proportion of pixels lost to an other class of a different crop/nocrop kind.")
    #    doc2.write("\\subsection{Results}")
    #    data1 = OptPerClassStatList["_OSORE"].NegMRMoy[Idx]
    #    data2 = RadPerClassStatList[""].NegMRMoy[Idx]
    #    data3 = FusPerClassStatList["_Fusion_OSORE"].NegMRMoy[Idx]
    #    PerClassTableau(cname,data1,data2,data3,DoY, title +": Negative miss rate averaged over 10 runs","Min")
    #
    #    # Precision per class
    #    doc2.write("\\section{Pecision per class}")
    #    doc2.write("\\subsection{Definition}")
    #    doc2.write("Precision per class: Proportion of pixels correctly classify among the predicted label.")
    #    doc2.write("\\subsection{Results}")
    #    data1 = OptPerClassStatList["_OSORE"].PreMoy[Idx]
    #    data2 = RadPerClassStatList[""].PreMoy[Idx]
    #    data3 = FusPerClassStatList["_Fusion_OSORE"].PreMoy[Idx]
    #    PerClassTableau(cname,data1,data2,data3,DoY, title +": Precision per class averaged over 10 runs","Max")
    # 
    #    # Positive False detection rate
    #    doc2.write("\\newpage")
    #    doc2.write("\\section{Positive false detection rate}")
    #    doc2.write("\\subsection{Definition}")
    #    doc2.write("Positive false detection rate: Proportion of pixels incorrectly classify to an other class that is however of the same crop/nocrop label.")
    #    doc2.write("\\subsection{Results}")
    #    data1 = OptPerClassStatList["_OSORE"].PosFDRMoy[Idx]
    #    data2 = RadPerClassStatList[""].PosFDRMoy[Idx]
    #    data3 = FusPerClassStatList["_Fusion_OSORE"].PosFDRMoy[Idx]
    #    PerClassTableau(cname,data1,data2,data3,DoY, title +": Positive false detection rate averaged over 10 runs","Min")
    #
    #    # Ngative miss rate per class
    #    doc2.write("\\newpage")
    #    doc2.write("\\section{Negative false detection rate}")
    #    doc2.write("\\subsection{Definition}")
    #    doc2.write("Negative false detection rate: Proportion of pixels incorrectly classify to an other class that is of a different crop/nocrop label.")
    #    doc2.write("\\subsection{Results}")
    #    data1 = OptPerClassStatList["_OSORE"].NegFDRMoy[Idx]
    #    data2 = RadPerClassStatList[""].NegFDRMoy[Idx]
    #    data3 = FusPerClassStatList["_Fusion_OSORE"].NegFDRMoy[Idx]
    #    PerClassTableau(cname,data1,data2,data3,DoY, title +": Negative fasle detection rate averaged over 10 runs","Min")
    # 
    #    doc2.close()
     
    
    # Construction of the Fscore Plots
    # Loop over classes
    ci = 0 
    cf = NbClasses

    for c in range(ci,cf):
        print "### %s ###"%(classname[int(OptFscore[0][c][0])])
        fig = plt.figure(figsize=(12,6))
        plt.grid()
        #plt.ylim(0, 1)
        plt.title(PlotTitle + classname[int(OptFscore[0][c][0])])
        cname = classname[int(OptFscore[0][c][0])].replace(" ", "_")
        namepng="fig/%s%s-Fscore-%s.png"%(title,TypeName[croptype],cname)
        plt.xlabel('DoY', fontsize = 13)
        plt.ylabel('Average F-Score', fontsize = 13)
        i=0 # Color index #
        if(OpticalTest):
            print "    *** Optical ***"
            StatList = []
            date = []
            for btype in OptBtypeList:
                print "    ",bandstype[btype]
                OptStat = OptStatList[btype]
                fscore = OptStat.FSmoy
                fscoreOPT = fscore
                OptOA = OptStat.OAmoy
                OptInt = OptStat.OAint
                fsbar = OptStat.FSint
                Optdatesidx = OptStat.dates   
                StatList.append(OptStat)
                date.append(OptDate[Optdatesidx])
                plt.title(PlotTitle + classname[int(OptFscore[0][c][0])])
                plt.errorbar(OptDate[Optdatesidx],fscore[:,c,1],fsbar[:,c,1],label= bandstype[btype],color = col[i])
                plt.legend(bbox_to_anchor=(xl, yl),bbox_transform=plt.gcf().transFigure)
                i += 1
            sorteddates = OptDate[Optdatesidx]
            DoY = [DoYDic[k] for k in sorteddates]
            plt.xticks(sorteddates,DoY,rotation = 90,fontsize = 10)
        if(RadarTest):
            print "    *** Radar ***"
            #plt.xticks(RadDate[::2], RadDoYlabel[::2],rotation = 90,fontsize = 10)
            for btype in RadBtypeList:
                print "    ",bandstype[btype]
                RadStat = RadStatList[btype]
                fscore = RadStat.FSmoy
                fscoreRAD = fscore
                RadOA = RadStat.OAmoy
                RadInt = RadStat.OAint
                fsbar = RadStat.FSint
                Raddatesidx = RadStat.dates
                StatList.append(RadStat)
                date.append(RadDate[Raddatesidx])
                plt.errorbar(RadDate[Raddatesidx],fscore[:,c,1],fsbar[:,c,1],label="Radar",color = col[i])
                plt.legend(bbox_to_anchor=(xl, yl),bbox_transform=plt.gcf().transFigure)
                i += 1
        if(FusionTest):
            print "    *** Fusion ***"
            for btype in FusBtypeList:
                print "    ",bandstype[btype]
                FusStat = FusStatList[btype]
                fscore = FusStat.FSmoy
                fscoreFUS = fscore
                FusOA = FusStat.OAmoy
                FusInt = FusStat.OAint
                fsbar = FusStat.FSint
                Fusdatesidx = FusStat.dates
                StatList.append(FusStat)
                
                # ATTENTION CHOIX DE DATE DE REFERENCE
                # Ref Optical #
                FusDate = OptDate
                date.append(OptDate[Fusdatesidx])
                plt.title(PlotTitle + classname[int(OptFscore[0][c][0])])
                plt.errorbar(FusDate[Fusdatesidx],fscore[:,c,1],fsbar[:,c,1],label= bandstype[btype],color = col[i])
                plt.legend(bbox_to_anchor=(xl, yl),bbox_transform=plt.gcf().transFigure)
                i += 1
            #plt.xticks(OptDate[datesidx], [OptDoYlabel[ii] for ii in datesidx],rotation = 90,fontsize = 10)
    
        if(png == 0): pdf.savefig(fig,bbox_inches='tight')
        else:plt.savefig(namepng)
    if(png == 0):pdf.close()
    
    
    Idx = [2,4,6]
    print RadDoYlabel
    for i in Idx:
        print "RAD:",RadDoYlabel[Raddatesidx[i]]
        print "OPT:",OptDoYlabel[Optdatesidx[i]]
    
    
    
    print StatList[3].BinRECmoy[:,0][[2,4,6]],StatList[3].BinPREmoy[:,0][[2,4,6]]
    print StatList[3].BinRECint[:,0][[2,4,6]],StatList[3].BinPREint[:,0][[2,4,6]]
    print StatList[0].BinRECmoy[:,0][[2,4,6]],StatList[0].BinPREmoy[:,0][[2,4,6]]
    print StatList[0].BinRECint[:,0][[2,4,6]],StatList[0].BinPREint[:,0][[2,4,6]]
    print "***********"
    print StatList[3].BinRECmoy[:,1][[2,4,6]],StatList[3].BinPREmoy[:,1][[2,4,6]]
    print StatList[3].BinRECint[:,1][[2,4,6]],StatList[3].BinPREint[:,1][[2,4,6]]
    print StatList[0].BinRECmoy[:,1][[2,4,6]],StatList[0].BinPREmoy[:,1][[2,4,6]]
    print StatList[0].BinRECint[:,1][[2,4,6]],StatList[0].BinPREint[:,1][[2,4,6]]
    
    
    
    # Construction of the crop mask statistic #
    if(croptype == 0):
        xl = 0.99 
        yl = 0.25
        labelname = ["S2AGRI","Opt + OSORE", "Radar","Fusion"]
     
        #labelname = ["Optical", "Opt + OSORE", "Radar"]
        #labelname = ["Opt", "Radar","Fusion"]
        #labelname = ["Opt", "Opt + OSORE", "Radar", "Fusion","Fusion OSORE"]
        print "Plot CropMask Statistic"
        CropMaskPDF = outputdir+"Plot-PrRe-CropMask.pdf"
        if(png == 0): pdf = matplotlib.backends.backend_pdf.PdfPages(CropMaskPDF) 
        namepng="fig/%s%s-Pr0.png"%(title,TypeName[croptype])
        #plt.ylim(0, 1) 
        print namepng
        fig = plt.figure(figsize=(10,6))
        plt.grid()
        for i,stat in enumerate(StatList):
            plt.title(title + ": No Crop Precision")
            plt.xlabel('doy', fontsize = 13)
            plt.ylabel('NoCropPre.', fontsize = 13)
            plt.errorbar(date[i],stat.BinPREmoy[:,0],stat.BinPREint[:,0],label = labelname[i],color = col[i])
            plt.legend(bbox_to_anchor=(xl, yl))
            plt.xticks(sorteddates,DoY,rotation = 90,fontsize = 10)
        if(png == 0): pdf.savefig(fig,bbox_inches='tight')
        else:plt.savefig(namepng,bbox_inches='tight')
        
        namepng="fig/%s%s-Re0.png"%(title,TypeName[croptype])
        fig = plt.figure(figsize=(10,6))
        plt.grid()
        for i,stat in enumerate(StatList):
            plt.title(title + ": No Crop Recall")
            plt.xlabel('doy', fontsize = 13)
            plt.ylabel('NoCropRec.', fontsize = 13)
            #print "Re0 = ", stat.BinRECmoy[:,0]
            plt.errorbar(date[i],stat.BinRECmoy[:,0],stat.BinRECint[:,0],label = labelname[i],color = col[i])
            plt.legend(bbox_to_anchor=(xl, yl))
            plt.xticks(sorteddates,DoY,rotation = 90,fontsize = 10)
        if(png == 0): pdf.savefig(fig,bbox_inches='tight')
        else:plt.savefig(namepng,bbox_inches='tight')
      
        namepng="fig/%s%s-Pr1.png"%(title,TypeName[croptype])
        fig = plt.figure(figsize=(10,6))
        plt.grid()
        for i,stat in enumerate(StatList):
            plt.title(title + ": Crop Precision")
            plt.xlabel('doy', fontsize = 13)
            plt.ylabel('CropPre.', fontsize = 13)
            plt.errorbar(date[i],stat.BinPREmoy[:,1],stat.BinPREint[:,1],label = labelname[i],color = col[i])
            plt.legend(bbox_to_anchor=(xl, yl))
            plt.xticks(sorteddates,DoY,rotation = 90,fontsize = 10)
        if(png == 0): pdf.savefig(fig,bbox_inches='tight')
        else:plt.savefig(namepng,bbox_inches='tight')
      
        namepng="fig/%s%s-Re1.png"%(title,TypeName[croptype])
        fig = plt.figure(figsize=(10,6))
        plt.grid()
        for i,stat in enumerate(StatList):
            plt.title(title + ": Crop Recall")
            plt.xlabel('doy', fontsize = 13)
            plt.ylabel('CropRec.', fontsize = 13)
            #print "Re1 = ", stat.BinRECmoy[:,1]
            plt.errorbar(date[i],stat.BinRECmoy[:,1],stat.BinRECint[:,1],label = labelname[i],color = col[i])
            plt.legend(bbox_to_anchor=(xl, yl))
            plt.xticks(sorteddates,DoY,rotation = 90,fontsize = 10)
        if(png == 0): pdf.savefig(fig,bbox_inches='tight')
        else:plt.savefig(namepng,bbox_inches='tight')
     
        #plt.errorbar(OptDate[Optdatesidx],OptStat.BinRECmoy[:,0],OptStat.BinRECint[:,0],label= "NoCrop Rec.",color = col[0],)
        # plt.errorbar(OptDate[Optdatesidx],OptStat.BinPREmoy[:,1],OptStat.BinPREint[:,1],label= "Crop Pre.",color = col[1])
        #plt.errorbar(OptDate[Optdatesidx],OptStat.BinRECmoy[:,1],OptStat.BinRECint[:,1],label= "Crop Rec",color = col[1],)
        #plt.legend(bbox_to_anchor=(1.3, 0.5))
        if(png == 0): pdf.close()
    
    # construction of the Confusion Plots#
    # Set dates indices
    if(title=="France2016"):
        Idx = [2,5,8]
    else:
        Idx = [2,4,6]
    
    
    width = 50 
    print "Plot Confusions"
    pdf = matplotlib.backends.backend_pdf.PdfPages(ConfusionsPDF)
    plt.ylim([0,100])
    
    for cj in range(FusStat.Conf.shape[1]):
    #for cj in range(1):
      fig = plt.figure(figsize=(12,6))
      conftitle = "%s: False Positives for %s"%(title,classname[int(OptFscore[0][cj][0])])
      print conftitle
      plt.title(conftitle)
      #namepng="fig/%s%s-Fscore-%s.png"%(title,TypeName[croptype],classname[int(OptFscore[0][c][0])])
      plt.ylim(0, 100)
      plt.xlabel('DoY', fontsize = 13)
      plt.ylabel('Proportion', fontsize = 13)
      plt.xticks(sorteddates,DoY,rotation = 90,fontsize = 10)
      for d in range(FusStat.Conf.shape[0]):
        ref = 0
        for ci in range(FusStat.Conf.shape[2]):
          ClassName = classname[int(OptFscore[0][ci][0])]
          val = 100.0*FusStat.Conf[d,ci,cj]/np.sum(FusStat.Conf[d,:,cj])
          if(d == 0):
            plt.bar(FusDate[Fusdatesidx[d]],val, width, color = colordic[int(OptFscore[0][ci][0])], label = ClassName,bottom = ref)
            plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0)
          else: plt.bar(FusDate[Fusdatesidx[d]],val, width, color = colordic[int(OptFscore[0][ci][0])],bottom = ref)
          ref = ref + val
      pdf.savefig(fig,bbox_inches='tight')
      
      fig = plt.figure(figsize=(12,6))
      conftitle = "%s: False Negatives for %s"%(title,classname[int(OptFscore[0][cj][0])])
      print conftitle
      plt.title(conftitle)
      #namepng="fig/%s%s-Fscore-%s.png"%(title,TypeName[croptype],classname[int(OptFscore[0][c][0])])
      plt.ylim(0, 100)
      plt.xlabel('DoY', fontsize = 13)
      plt.ylabel('Proportion', fontsize = 13)
      plt.xticks(sorteddates,DoY,rotation = 90,fontsize = 10)
      for d in range(FusStat.Conf.shape[0]):
        ref = 0
        for ci in range(FusStat.Conf.shape[2]):
          ClassName = classname[int(OptFscore[0][ci][0])]
          val = 100.0*FusStat.Conf[d,cj,ci]/np.sum(FusStat.Conf[d,cj,:])
          if(d == 0):
            plt.bar(FusDate[Fusdatesidx[d]],val, width, color = colordic[int(OptFscore[0][ci][0])], label = ClassName,bottom = ref)
            plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0)
          else: plt.bar(FusDate[Fusdatesidx[d]],val, width, color = colordic[int(OptFscore[0][ci][0])],bottom = ref)
          ref = ref + val
      pdf.savefig(fig,bbox_inches='tight')
    pdf.close()
    
    # Construction of the Fscore Tables
    # Part that should be improve because so far "the strategic" dates have been fixed "by hand"
    doc = latex.latexdoc(TablePDF)
    doc.title(title)
    doc.write("\\section{Averaged F-score}")
    DoY = [RadDoYlabel[Raddatesidx[Idx[0]]],RadDoYlabel[Raddatesidx[Idx[1]]],RadDoYlabel[Raddatesidx[Idx[2]]]]
    printFscore(fscoreRAD,RadOA,RadInt,fscoreOPT,OptOA,OptInt,fscoreFUS,FusOA,FusInt,Idx,Idx,DoY,classname,doc,"Averaged F-score over 10 runs.")
    doc.close()
    
    
    ###############################################
    # Calculate statistic after classes gathering #
    ###############################################
    
    OptStatList = {}
    OptPerClassStatList = {}
    RadStatList = {}
    RadPerClassStatList = {}
    FusStatList = {}
    FusPerClassStatList = {}
    for btype in OptBtypeList:
        print "Get Optical Measures"
        OptStatList[btype] = FscoreResultsGather(OptDirectory,RFDirOpt,NbOptDates,Nbtirages,btype,cropmask,imperfectCM,gatheringSet)
        #if(croptype == 0): OptPerClassStatList[btype] = PerClassResults(OptDirectory,RFDirOpt,NbOptDates,Nbtirages,btype,cropmask)
    
    for btype in RadBtypeList:
        print "Get Radar Measures"
        RadStatList[btype] = FscoreResultsGather(RadDirectory,RFDirRad,NbRadDates,Nbtirages,btype,cropmask,imperfectCM,gatheringSet)
        #if(croptype == 0): RadPerClassStatList[btype] = PerClassResults(RadDirectory,RFDirRad,NbRadDates,Nbtirages,btype,cropmask)
    
    for btype in FusBtypeList:
        print "Get Fusion Measures"
        FusStatList[btype] = FscoreResultsGather(FusDirectory,RFDirOpt,NbOptDates,Nbtirages,btype,cropmask,imperfectCM,gatheringSet)
        #if(croptype == 0): FusPerClassStatList[btype] = PerClassResults(FusDirectory,RFDirOpt,NbOptDates,Nbtirages,btype,cropmask)
    
    # FRANCE
    classname[1] = "Winter Crops"
    classname[2] = "Summer Crops"
    
    # ITALY
    #classname[1] = "Winter Crops"
    #classname[2] = "Leguminous Crops"
    #classname[3] = "Vegetables"
    
    # SPAIN
    #classname[1] = "Winter Crops"
    #classname[2] = "Summer Crops"
    #classname[3] = "Leguminous Crops"
    #classname[4] = "Vegetables"
    #classname[5] = "Roots"
    
    
    wsc = ""
    if(title=="Spain2017"): wsc = "_WSC"
    
    print RadDoYlabel
    for i in Idx:
        print "RAD:",RadDoYlabel[Raddatesidx[i]]
        print "OPT:",OptDoYlabel[Optdatesidx[i]]
    
    
    doc = latex.latexdoc(TablePDFgather)
    doc.title(title)
    doc.write("\\section{Averaged F-score}")
    DoY = [RadDoYlabel[Raddatesidx[Idx[0]]],RadDoYlabel[Raddatesidx[Idx[1]]],RadDoYlabel[Raddatesidx[Idx[2]]]]
    printFscoreNew(RadStatList["_CROP%s"%(wsc)],OptStatList["_CROP_OSORE%s"%(wsc)],FusStatList["_CROP_Fusion_OSORE%s"%(wsc)],Idx,Idx,DoY,classname,doc,"Averaged F-score over 10 runs.")
    doc.close()
    
    # Si probleme "dictionary", faire la boucle sur c en plusieurs fois
    # puis concatene fichiers pdf avec pdfunite (pdfunite 1.pdf 2.pdf out.pdf)
