#!/usr/bin/env python
#-*- coding: utf-8 -*-
######################################
#          CESBIO 2017-2019          #
#       Creation: Ludo 06/06/2019    #
######################################

import os
import numpy as np
import ConfigParser
import json
from dates import GetDates,GetDatesSentinel
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as mpdf
import datetime
import PlotStatisticsNew as psn
import PlotConfusionMatrix as pcm
#import PlotVariablesImportance as pvi
import tolatex.tolatex as latex

def printError():
    print("\033[0;31mERROR:")
    os.system("tput init") # Reset color

def Export(fileconf):
    """ Export Statistical Results """
    print "Exporting results plots and tables..."
    config = ConfigParser.ConfigParser()
    namecfg = config.read(fileconf)

    try:
        title = config.get("Global","Title") 
        NbRun = int(config.get("Global","NbRun"))
        outputdir = config.get("Global","OutputDirectory")
    
        fname = config.get("Global","Name")
      
        SectionsList = config.sections()
        SectionsList.remove("Global")
    except:
        printError()
	print "It seems that the results configuration file does not exist or that its structure is incorrect. Please check."
	quit()

    
    CropTypePDF   = outputdir + title + "-Plot-RePrCropMask-%s.pdf"%(fname)
    FscorePDF     = outputdir + title + "-Plot-Fscore-AllClasses-%s.pdf"%(fname)
    FscoreBAR     = outputdir + title + "-Bar-Fscore-AllClasses-%s.pdf"%(fname)
    TablePDF      = outputdir + title + "-Table-Fscore-AllClasses-%s"%(fname)
    ConfusionsPDF = outputdir + title + "-Plot-Confusions-AllClasses-%s-%s.pdf"%(fname,"%s")
    VIPDF = outputdir + title + "-Plot-VariablesImportance-AllClasses-%s-%s.pdf"%(fname,"%s")

    PrReCM = mpdf.PdfPages(CropTypePDF)
    pdf = mpdf.PdfPages(FscorePDF)
    bar = mpdf.PdfPages(FscoreBAR)
    vi = mpdf.PdfPages(VIPDF)

    ListBand = []
    ListName = []
    ListDates = []
    ListStyle = []
    ListStat = []
    ListValDates = []
    ListTableDates = []
    ListQTable = []

    for section in SectionsList:
      # TODO: Ajouter comme param
      croptype = 0 #TODO  
      imperfectCM = False
      TypeName=["","-CROP","CTAfterCM"]
      png = 0
      xl = 0.9 
      yl = 0.3

      # Load config parameters
      workdir = config.get(section,"WorkDirectory")
      datadir = config.get(section,"WorkDirectory")+"/AllTiles/"+config.get(section,"DataDirectory")+"/"
      RFDir = config.get(section,"ClassifDirectory")
      band = config.get(section,"Label")
      name = config.get(section,"Name")
      DatesFile = workdir + "/WorkFiles/" + config.get(section,"DatesFile")
      ValDates = json.loads(config.get(section,"ValidationDates"))
      style = json.loads(config.get(section,"Style"))
      QTable = config.getboolean(section,"Table")
      if(QTable):
        try:
          TableDates = json.loads(config.get(section,"TableDates"))
          ListTableDates.append(TableDates)
        except:
          print "TableDates parameter need to be provided in section %s"%(section)
          quit()

      try:
	      Qnorm = config.getboolean(section,"Normalized")
      except:
	      Qnorm = False


      # Path and directory settings
      Dates,DoYlabel = GetDates(DatesFile)
      NbDates = len(Dates)
      Dates = np.array(Dates)
      DoYDic = {}
      DoYList = []
      for i,d in enumerate(Dates):
          DoYDic[d]=DoYlabel[i]
          DoYList.append(DoYlabel[i])
       
      #classname,gatheringSet,colordic,cropmask = ConstrucDics(workdir,croptype)
      classname,colordic,cropmask = ConstrucDics(workdir,croptype)


      # TODO          
      # Set variables for "Crop type after crop mask"
      #TypeName=["","-CROP","CTAfterCM"]
      #if (png==1):os.system("mkdir -p fig")

      PlotTitle =  title + ": All Classes Fscore: "

      # Get all statistic per classification type 
      Stat = psn.FscoreResults(datadir,RFDir,NbDates,NbRun,band,cropmask,imperfectCM,Qnorm)

      # Fill param list per configuiration section
      ListBand.append(band)
      ListName.append(name)
      ListStyle.append(style)
      ListValDates.append(ValDates)
      ListDates.append([NbDates,Dates,DoYDic,DoYList,TableDates])
      ListStat.append(Stat)
      ListQTable.append(QTable)


      # Get out Number of classe
      Fscore = Stat.FSmoy
      NbClasses= len(Fscore[0])
 
       
      # Construction of the Fscore Plots
      # Loop over classes
    
    ND = len(ListTableDates[0])
    NS = len(ListQTable)
    tdata=[]
    bardata=[]
    errdata=[]
    ci = 0 
    cf = NbClasses
    ClassNameList = []

    print "*** Export Average Fscore evolution ***"
    for c in range(ci,cf):
        print "   - %s"%(classname[int(Fscore[0][c][0])])
        ClassNameList.append(classname[int(Fscore[0][c][0])])
        figPlot = plt.figure(figsize=(12,6))

        axPlot = figPlot.add_subplot(111)

        axPlot.grid()
        #plt.ylim(0, 1)

        axPlot.set_title(PlotTitle + classname[int(Fscore[0][c][0])])
        cname = classname[int(Fscore[0][c][0])].replace(" ", "_")
        namepng="fig/%s%s-Fscore-%s.png"%(title,TypeName[croptype],cname)
        axPlot.set_xlabel('DoY', fontsize = 13)
        axPlot.set_ylabel('Average F-Score (%)', fontsize = 13)


    # ListBand = []
    # ListName = []
    # ListDates = []
    # ListCol = []
    # ListStat = []
        ListFS = []
        ListErr = []
        line = []
        for s,section in enumerate(SectionsList):
          #print "   *** ",ListName[s]," ***"
          #print "   *** ",ListBand[s]," ***"
          # Get Statistical variable for all sections 
          Stat = ListStat[s]
          fscore = Stat.FSmoy
          OA = Stat.OAmoy
          Int = Stat.OAint
          fsbar = Stat.FSint
          datesidx = Stat.dates   
          
          # Plot variables per sections
          axPlot.set_title(PlotTitle + classname[int(Fscore[0][c][0])])
          axPlot.errorbar(ListDates[s][1][datesidx],100.0*fscore[:,c,1],100.0*fsbar[:,c,1],capsize=2,lw=0.5,label = ListName[s], color = ListStyle[s][0], ls = ListStyle[s][1])
          #axPlot.legend(bbox_to_anchor=(xl, yl),bbox_transform=plt.gcf().transFigure)
          axPlot.legend(loc=0)
          sorteddates = ListDates[s][1][datesidx]
          DoY = [ListDates[s][2][k] for k in sorteddates]
          axPlot.set_xticks(sorteddates)
          axPlot.set_xticklabels(DoY,rotation = 90,fontsize = 10)

          # Get fscore per section at the chosend date to construct table
          TableDatesIdx = [ListValDates[s].index(i) for i in ListTableDates[s]]
          fs = 100.0*fscore[TableDatesIdx,c,1]
          err = 100.0*fsbar[TableDatesIdx,c,1]
          if(ListQTable[s]):
		  ListFS.append(fs)
		  ListErr.append(err)

        #print ListFS
        # Transpose List of fscore from Section-Date to Date Section
        line = [classname[int(Fscore[0][c][0])]]
        barline = []
        errline = []
        for i in range(ND): 
          for Fs in ListFS:
            #line.append(Fs[i])
            line.append(Fs[i]/100.0)
            barline.append(Fs[i])
          for Err in ListErr:
            errline.append(Err[i])

        tdata.append(line)
        bardata.append(barline)
        errdata.append(errline)

        if(png == 0):
	  pdf.savefig(figPlot,bbox_inches='tight')
        #else:plt.savefig(namepng)
    if(png == 0):
	    pdf.close()   

     
########## Table ##############

    print "*** Export Average Fscore table ***"
    doc = latex.latexdoc(TablePDF)
    doc.title(title)

    # Construct table format 
    tf="|l|"
    for i in range(ND): 
      dc="|"
      for j in range(NS):
        dc=dc+"c|"
      tf=tf+dc
    
    # Construct band header
    head = ["Class"]
    multisize = [1]
    for i in range(ND):
      multisize.append(len(SectionsList))
      for s,section in enumerate(SectionsList):
        if(ListQTable[s]): head.append(section)
    #print "multisize",multisize

    
    # Construct Date header
    mc = [""]
    for d in ListDates[0][4]:
      mc.append(ListDates[0][3][d-1])

    #print tdata
    # Add bold for min/max
    mark = "Max"
    for i in range(len(tdata)):
      line = tdata[i]
      for d in range(ND):
	currentline = line[NS*d+1:NS*d+NS+1]
        if(mark=="Max"):idx = line.index(max(currentline))
        elif(mark=="Min"):idx = line.index(min(currentline))
	try:
		tdata[i][idx] =  "\\bf{%.2f}"%(line[idx])
	except:
		raise
    #tdata.insert(0,head)
    alldata = tdata

    cap = "%s: F-Score averaged over %d random runs"%(title,NbRun)
    doc.tableau([[head],alldata],tf, multicol = mc, size = multisize,caption = cap + " (%simum are indicated with bold font)."%(mark))
    doc.write("\\footnotetext{Documents created automatically by the SenSAgri Classification Chain on \\date{\\oldstylenums{\\today}}}")
    doc.close()


################# BAR #################"" 

    print "*** Export Average Fscore barplot ***"
    bardata = np.asarray(bardata)
    errdata = np.asarray(errdata)

    width = 1.0/float(NS+1)
    x = np.arange(NbClasses)
    dt = len(ListDates[0][4])
    for i,d in enumerate(ListDates[0][4]):
            figBar  = plt.figure(figsize=(12,6))
            axBar  = figBar.add_subplot(111)

            axBar.set_ylim(0, 110)
            axBar.grid()
            axBar.set_title(PlotTitle)
            axBar.set_xlabel('Class', fontsize = 13)
            axBar.set_ylabel('Average F-Score (%)', fontsize = 13)
	    axBar.set_title(PlotTitle + ListDates[0][3][d-1])
            for s,section in enumerate(SectionsList):
		    Stat = ListStat[s]
	            fscore = Stat.FSmoy
		    fsbar = Stat.FSint
		    
		    fscoreBar = bardata[:,s+i*NS] 
		    errorBar = errdata[:,s+i*NS]
		    shift = width*s
		    axBar.bar(x+shift,fscoreBar,width*0.9,yerr=errorBar,label = ListName[s],color = ListStyle[s][0], ecolor = "k")
          	    axBar.set_xticks(range(NbClasses))
           	    axBar.set_xticklabels(ClassNameList,rotation = 90,fontsize = 10)

            axBar.legend()
            bar.savefig(figBar,bbox_inches='tight')
   
    bar.close()   

################# Crop Mask ################# 

    print "*** Export Crop Mask measures evolution ***"
    measuresName = {"Re0":"NoCrop Recall","Pr0":"NoCrop Precision","Re1":"Crop Recall","Pr1":"Crop Precision"}
    measures     = {"Re0":"BinRECmoy","Pr0":"BinPREmoy","Re1":"BinRECmoy","Pr1":"BinPREmoy"}
    measuresInt  = {"Re0":"BinRECint","Pr0":"BinPREint","Re1":"BinRECint","Pr1":"BinPREint"}
    iscrop =  {"Re0":0,"Pr0":0,"Re1":1,"Pr1":1}

    meas = {}
    for meas in ["Re0","Pr0","Re1","Pr1"]:
	    figCM  = plt.figure(figsize=(12,6))
            axCM  = figCM.add_subplot(111)
            #axCM.set_ylim(0, 100)
            axCM.grid()
            axCM.set_title(PlotTitle)
            axCM.set_xlabel('DoY', fontsize = 13)
            axCM.set_ylabel('Average %s (%%)'%(measuresName[meas]), fontsize = 13)
            axCM.set_title(PlotTitle + measuresName[meas])
            for s,section in enumerate(SectionsList):
                Stat = ListStat[s]
                datesidx = Stat.dates
                ##Meas = Stat.BinRECmoy[:,0]
                Meas = getattr(Stat,measures[meas])[:,iscrop[meas]]
                IntMeas = getattr(Stat,measuresInt[meas])[:,iscrop[meas]]
                axCM.errorbar(ListDates[s][1][datesidx],100.0*Meas,100.0*IntMeas,capsize=2,lw=0.5,label = ListName[s], color = ListStyle[s][0], ls = ListStyle[s][1])
                axCM.legend()
                sorteddates = ListDates[s][1][datesidx]
                DoY = [ListDates[s][2][k] for k in sorteddates]
                axCM.set_xticks(sorteddates)
                axCM.set_xticklabels(DoY,rotation = 90,fontsize = 10)
                axCM.legend()
            PrReCM.savefig(figCM,bbox_inches='tight')
    PrReCM.close()   

#################### Confusion ##############################

# TODO

    print "*** Export Confusion Matrices ***" 
    for s,section in enumerate(SectionsList):
      print "  - %s"%(section)
      TableDatesIdx = [ListValDates[s].index(i) for i in ListTableDates[s]]
      oalist = 100.0*ListStat[s].OAmoy[TableDatesIdx]
      intlist = 100.0*ListStat[s].OAint[TableDatesIdx]
      reclist = 100.0*ListStat[s].RECmoy[TableDatesIdx,:,1]
      fslist = 100.0*ListStat[s].FSmoy[TableDatesIdx,:,1]
      prelist = 100.0*ListStat[s].PREmoy[TableDatesIdx,:,1]
      cmlist = ListStat[s].Conf[TableDatesIdx]
      datesList = [ListDates[s][3][ListTableDates[s][i]-1] for i,fs in enumerate(fslist)]

      
      pcm.PlotCM(cmlist,oalist,intlist,reclist,prelist,fslist,datesList,ClassNameList,ConfusionsPDF,title,section)


def ConstrucDics(workdir,croptype):
   """ Construct Output Directories  """  
   classfile = workdir + "/WorkFiles/Classes.csv" 
   #gatherfile = workdir + "/WorkFiles/Classes_Gathered.csv" 
   binaryfile = workdir + "/WorkFiles/Binary.txt" 
   colorfile = workdir + "/WorkFiles/Legend.csv" 
   #outputdir = workdir + "/Statistics/" 
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
      
  # # Construct Class gathering dictionary 
  # with open(gatherfile, "r") as ins:
  #     gatheringSet = {}
  #     for line in ins:
  #         couple = line.split(":")
  #         gatheringSet[int(couple[0])] = int(couple[1])
  #     gatheringSet[0] = 0
     
   # Construct Color dictionary
   #with open(colorfile, "r") as ins:
   #    colordic = {0:"#000000"}
   #    for line in ins:
   #        couple = line.split(":")
   #        colordic[int(couple[0])] = couple[1][:-1] # [:-1] to remove the "\n" at the end of the line
   # # Fix cropmask "mask"

   colordic = {0:0}

   if croptype == 0:
       cropmask = np.loadtxt(binaryfile)
       # filtered
       cropmask = cropmask[:-1].copy()
   else: cropmask = np.array([1])
   #return classname,gatheringSet,colordic,cropmask
   return classname,colordic,cropmask
