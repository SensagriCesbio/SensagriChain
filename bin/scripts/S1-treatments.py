#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import sys
import csv

tiles = ["30TUM","30TUL"]
ConIntImage = {"30TUM":"S1_T30TUM_2018_Concatenate.tif","30TUL":"S1_T30TUL_2018_Concatenate.tif"}
path = "/tmpdir/larnaud/2018/Spain2018/WorkFiles/Images/S1-"
srun = "srun --hint=nomultithread -N1 -n1 -c1 --exclusive "
header = """#!/bin/bash
#SBATCH --job-name=S1-Stacks
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --ntasks-per-node=3
#SBATCH --mem-per-cpu=60G
#SBATCH --time=2-12:00:00
#SBATCH --output="S1-Stacks-%j.out"
#SBATCH --error="S1-Stacks-%j.err"
#SBATCH --mail-user=user@server.com
#SBATCH --mail-type=ALL

# Log function #
printlog () {
  echo LOG $(date +'%d/%m/%Y-%H:%M:%S'): $1
}
printdebug () {
  echo DBG $(date +'%d/%m/%Y-%H:%M:%S'): $1
}
# Load MPI module #
module purge
module load intel/18.2 intelmpi/18.2 chdb/1.0

export OTB_MAX_RAM_HINT=45000
dirname=/tmpdir/larnaud/2018/France2018

"""

outfile = open("s1stacks.sh","w")
outfile.write(header)

QCalculateDates = False

Qprepare = True
Qconca = True
Qinter = False
Qsuper = True

if(Qprepare):
    print("*** Prepartion directory **")
    datesList = []
    imagesList = {}
    for t in tiles:
        print t
        os.system("mkdir -p T%s"%(t))
	print "%s%s"%(path,t)
        if(QCalculateDates):
            os.system("ls %s%s > T%s/ListImages.txt"%(path,t,t))
            os.system("sort -n -k 1.22,1.29 T%s/ListImages.txt > T%s/Sorted_ListImages.txt"%(t,t))
            os.system("""cut -c11-29 T%s/Sorted_ListImages.txt | grep "vh" | cut -c 12- > T%s/S1_T%s_Dates.txt"""%(t,t,t))
        listfile = "T%s/Sorted_ListImages.txt"%(t)
        images = []
        dates = []
        with open(listfile, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                images.append(row[0])
                dates.append(row[0][21:29])
        datesList.append(dates)
        imagesList[t] = images
        print "Images number = ",len(images)

    #for i,d in enumerate(datesList[0]):
    #    print datesList[0][i],datesList[1][i],datesList[2][i],datesList[0][i]==datesList[1][i]
   
    for i,t in enumerate(tiles):
      for j,u in enumerate(tiles):
        if (j>i):
          if(datesList[i]==datesList[j]):
            print t,u,"   -> Same dates"
          else:
            print t,u,"   -> Different dates. Dates interpolation needed"

if(Qconca):
    print("*** Concatenation **")
    outfile.write("""printlog "Concatenation"\n""")
    for t in tiles:
        command = srun + "otbcli_ConcatenateImages -il"
        for im in imagesList[t]:
            command = command + " %s%s/"%(path,t) + im
        command = command + " -out T%s/S1_T%s_2018_Concatenate.tif > conc_%s.txt &\n"%(t,t,t)
        outfile.write(command)
    outfile.write("wait\n\n")
   
if(Qinter):
    # INTERPOLATION
    print "*** Interpolation ***"
    outfile.write("""printlog "Interpolation"\n""")
    command = srun + "Interpolation T30TYP/S1_T30TYP_2018_Concatenate.tif T30TYP/S1_T30TYP_Dates.txt T30TYP/S1_T30TYP_2018_Interpol.tif T31TCH/S1_T31TCH_Dates.txt 2 > inter_30TYP.txt &\n"
    #print command
    outfile.write(command)
    outfile.write("wait\n\n")

if(Qsuper):
    # SUPERIMPOSE
    print "*** Superimpose ***"
    outfile.write("""printlog "Superimpose"\n""")
    for t in tiles:
      inr = "/tmpdir/larnaud/2018/Spain2018/WorkFiles/Images/Sentinel2_ST_REFL_GAP_T%s_2018.tif"%(t)
      inm = "T%s/%s"%(t,ConIntImage[t])
      out = "T%s/S1_T%s_Superimposed.tif"%(t,t)
      command = srun + "otbcli_Superimpose -inr %s -inm %s -out %s -interpolator nn > super_T%s.txt &\n"%(inr,inm,out,t)
      outfile.write(command)
outfile.write("wait\n\n")
 
outfile.write("""printlog "Job Completed"\n""")
outfile.close()
