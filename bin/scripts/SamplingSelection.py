#!/usr/bin/env python
import os
import sys
import csv
import numpy as np
import infolog

def ReadFile(f):
	lines = []
	with open(f, "r") as ins:
		for row in ins:
			lines.append(row)
	return lines

def insert(source_str, insert_str, pos):
	    return source_str[:pos]+insert_str+source_str[pos:]

def class2bin(crops,c):
	return int(c in crops)

if (len(sys.argv) != 7):
	print("Usage",sys.argv[0] + " cropslistfile labels1 labels2 profile1 profiles2 outtId")
	quit()

if __name__ == '__main__':
  """
  Select sample between two set. Used for crop mix classification
  """

  log = infolog.infolog()
  log.msg("Loading inputs files")

  # Load input files
  crops = np.loadtxt(sys.argv[1])
  labels0 = np.loadtxt(sys.argv[2],dtype=np.int)
  labels1 = np.loadtxt(sys.argv[3],dtype=np.int)
  # Do not use numpy cause format change in the file between int and float
  profiles0 = ReadFile(sys.argv[4]) 
  profiles1 = ReadFile(sys.argv[5])
 
  log.msg("Selection")
  FileOutLabels = insert(sys.argv[3], "_" + sys.argv[6], -4)
  FileOutProfiles = insert(sys.argv[5], "_" + sys.argv[6], -4)

  # Selection 
  binary0 = np.array([class2bin(crops,c) for c in labels0])
  binary1 = np.array([class2bin(crops,c) for c in labels1])
  idx0 = np.where(binary0==0)
  idx1 = np.where(binary1==1)
  binaryOut = np.concatenate((binary0[idx0],binary1[idx1]))
  labelsOut = np.concatenate((labels0[idx0],labels1[idx1]))

  #print len(labels0),len(labels1),len(labels0[idx0]),len(labels1[idx1]),len(labelsOut)

  profilesOut = open(FileOutProfiles,"w")
  for i,row in enumerate(profiles0):
	  if(binary0[i] == 0):
		  profilesOut.write(row)

  for i,row in enumerate(profiles1):
	  if(binary1[i] == 1):
		  profilesOut.write(row)
  
  profilesOut.close()

  log.msg("Export to ouput files")
  np.savetxt(FileOutLabels,labelsOut, fmt='%d')
  #np.savetxt(FileOutProfiles,profilesOut, fmt='%d' ,delimiter=',')

  log.msg("Done")
