#!/usr/bin/python
# -*- coding: utf-8 -*- 
import os
import sys
import numpy
from numpy import loadtxt,savetxt
from dates import GetDates

sys.path.insert(0, 'ShapeFilesUtils')

#This function return the merged DoY vector containing the optical and radar acquisition dates
#  -- OptRadarDoY : It contains the DoY
#  -- idOptRadDoY : It indicates the sensor type : if the DoY index is a radar ( 1), otherwise optical (2) 

def getMergedDoYnum(RadarDoY, OpticalDoY):

	NbOpticalDates = len (OpticalDoY)
	NbRadarDates = len (RadarDoY)

	OptRadarDoY = numpy.zeros(NbOpticalDates + NbRadarDates)
	idOptRadDoY = numpy.zeros(NbOpticalDates + NbRadarDates)

	idOptical = 0
	idRadar = 0
	idOptRad = 0

	while (idOptical < NbOpticalDates and idRadar< NbRadarDates):

		if  RadarDoY[idRadar]< OpticalDoY [idOptical] :
			OptRadarDoY[idOptRad] = RadarDoY[idRadar]
			idOptRadDoY [idOptRad] = 1 
			idRadar = idRadar  + 1
		else:
			OptRadarDoY[idOptRad] = OpticalDoY [idOptical]
			idOptRadDoY[idOptRad] = 2
			idOptical = idOptical  + 1

		idOptRad = idOptRad  + 1

	if idOptical == NbOpticalDates:

		while (idOptRad < (NbOpticalDates + NbRadarDates)):
			OptRadarDoY[idOptRad] = RadarDoY[idRadar]
			idOptRadDoY [idOptRad] = 1 
			idRadar = idRadar  + 1
			idOptRad = idOptRad  + 1
	else:

		while (idOptRad < (NbOpticalDates + NbRadarDates)):
			OptRadarDoY[idOptRad] = OpticalDoY [idOptical]
			idOptRadDoY [idOptRad] = 2
			idOptical = idOptical  + 1
			idOptRad = idOptRad  + 1
	return OptRadarDoY, idOptRadDoY 

def getMergedDoY ( RadarDoYFile , OpticalDoYFile):

	RadarDoY = loadtxt(RadarDoYFile, delimiter=',')
	OpticalDoY = loadtxt(OpticalDoYFile, delimiter=',')

	NbOpticalDates = len (OpticalDoY)
	NbRadarDates = len (RadarDoY)

	OptRadarDoY = numpy.zeros(NbOpticalDates + NbRadarDates)
	idOptRadDoY = numpy.zeros(NbOpticalDates + NbRadarDates)

	idOptical = 0
	idRadar = 0
	idOptRad = 0

	while (idOptical < NbOpticalDates and idRadar< NbRadarDates):

		if  RadarDoY[idRadar]< OpticalDoY [idOptical] :
			OptRadarDoY[idOptRad] = RadarDoY[idRadar]
			idOptRadDoY [idOptRad] = 1 
			idRadar = idRadar  + 1
		else:
			OptRadarDoY[idOptRad] = OpticalDoY [idOptical]
			idOptRadDoY[idOptRad] = 2
			idOptical = idOptical  + 1

		idOptRad = idOptRad  + 1

	if idOptical == NbOpticalDates:

		while (idOptRad < (NbOpticalDates + NbRadarDates)):
			OptRadarDoY[idOptRad] = RadarDoY[idRadar]
			idOptRadDoY [idOptRad] = 1 
			idRadar = idRadar  + 1
			idOptRad = idOptRad  + 1
	else:

		while (idOptRad < (NbOpticalDates + NbRadarDates)):
			OptRadarDoY[idOptRad] = OpticalDoY [idOptical]
			idOptRadDoY [idOptRad] = 2
			idOptical = idOptical  + 1
			idOptRad = idOptRad  + 1
	return OptRadarDoY, idOptRadDoY

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------


def getSizeofData (idOptRadDoY, NbFeatRadar, NbFeatOpt):

	TF = 0
	for ids in idOptRadDoY:
		if  ids == 1 :
			TF = TF + NbFeatRadar
		else:
			TF = TF + NbFeatOpt
	return TF

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------


def getMergedData ( RadarDataFile , OpticalDataFile, idOptRadDoY, NbFeatRadar, NbFeatOpt):

	RadarData = loadtxt(RadarDataFile, delimiter=',')
	OpticalData = loadtxt(OpticalDataFile, delimiter=',')

	if (len(RadarData) != len(OpticalData)):
		raise Exception("Radar and optical do not have the same number of samples")

  	NbSamples = len(RadarData)
    	OutSet  = numpy.zeros(shape=(NbSamples , getSizeofData (idOptRadDoY, NbFeatRadar, NbFeatOpt)))
	
	for idSample in range(0,NbSamples):
		idCumDate = 0
		idRadar = 0
		idOpt = 0
		for idDate  in range(0, len(idOptRadDoY)):

			if  idOptRadDoY [idDate] == 1 :
				# Radar acquisition
				for idF in range(0,NbFeatRadar):
					OutSet[idSample][idCumDate + idF] = RadarData[idSample][idRadar*NbFeatRadar + idF]
				idRadar = idRadar  + 1
				idCumDate =idCumDate + NbFeatRadar

			else:
				# Optical acquisition
				for idF in range(0,NbFeatOpt):
					OutSet[idSample][idCumDate + idF] = OpticalData[idSample][idOpt*NbFeatOpt + idF]
				idOpt = idOpt  + 1
				idCumDate =idCumDate + NbFeatOpt

	return OutSet

def generate(home,RadarDoYFile,OpticalDoYFile,NbRadarFeatures,NbOpticalFeatures):
    
    OpticalDoY, OpticalDateTags = GetDates(home+OpticalDoYFile)
    RadarDoY, RadarDateTags = GetDates(home+RadarDoYFile)
    
    [OptRadDoY, idOptRadDoY] = getMergedDoYnum(RadarDoY, OpticalDoY)
    
    NbFeatures = numpy.empty_like(idOptRadDoY)
    if idOptRadDoY[0] == 1:
        NbFeatures[0] = NbRadarFeatures
    else:
        NbFeatures[0] = NbOpticalFeatures
    
    for i in range(1,len(idOptRadDoY)):
        if idOptRadDoY[i] == 1:
            NbFeatures[i] = NbFeatures[i-1] + NbRadarFeatures
        else: 
            NbFeatures[i] = NbFeatures[i-1] + NbOpticalFeatures
    # print idOptRadDoY
    # print OptRadDoY
    # print NbFeatures
    # print len(idOptRadDoY)
    # print len(OptRadDoY)
    savetxt(home + "/S1S2_DatesOrder.txt",idOptRadDoY, fmt='%d')
    savetxt(home + "/S1S2_Dates.txt",OptRadDoY, fmt='%d')
    savetxt(home + "/S1S2_NbFeaturesPerDates.txt",NbFeatures, fmt='%d')


if __name__ == "__main__":

        if len(sys.argv) != 5:
            print "Usage:"
            print "RadarDatesFile OpticalDatesFile NbrRadarFeatures NbrOpticalFeatures"
            quit()
	RadarDoYFile = sys.argv[1]
	OpticalDoYFile = sys.argv[2]
	NbRadarFeatures = int(sys.argv[3])
	NbOpticalFeatures = int(sys.argv[4])

        OpticalDoY, OpticalDateTags = GetDates(OpticalDoYFile)
        RadarDoY, RadarDateTags = GetDates(RadarDoYFile)

	[OptRadDoY, idOptRadDoY] = getMergedDoYnum(RadarDoY, OpticalDoY)


        NbFeatures = numpy.empty_like(idOptRadDoY)
        if idOptRadDoY[0] == 1:
            NbFeatures[0] = NbRadarFeatures
        else:
            NbFeatures[0] = NbOpticalFeatures

        for i in range(1,len(idOptRadDoY)):
            if idOptRadDoY[i] == 1:
                NbFeatures[i] = NbFeatures[i-1] + NbRadarFeatures
            else: 
                NbFeatures[i] = NbFeatures[i-1] + NbOpticalFeatures

        RadarDates = loadtxt(RadarDoYFile, delimiter=',')
	OpticalDates = loadtxt(OpticalDoYFile, delimiter=',')

        optidx = 0
        radidx = 0
        OptRadDates = []
        for idx in idOptRadDoY:
          if idx == 1:
            OptRadDates.append(RadarDates[radidx])
            radidx += 1
          elif idx ==2:
            OptRadDates.append(OpticalDates[optidx])
            optidx +=1

	print idOptRadDoY
#	print OptRadDoY
#	print NbFeatures
#        print len(idOptRadDoY)
#        print len(OptRadDoY)
        print OptRadDates

         
	
        savetxt("S1S2_DatesOrder.txt",idOptRadDoY, fmt='%d')
        savetxt("S1S2_Dates.txt",OptRadDoY, fmt='%d')
        savetxt("S1S2_DatesFormat.txt",OptRadDates, fmt='%d')
        savetxt("S1S2_NbFeaturesPerDates.txt",NbFeatures, fmt='%d')

