#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CesBIO 2017 #
import sys;
import os, osr
import shutil
import math
import random
from osgeo import gdal, ogr, osr
import time

def GetAreaOfPolygons(dataSource,field,cl):
    
    layer = dataSource.GetLayer()
    layer.SetAttributeFilter(field+" = "+str(cl))
    LayerArea = 0
    for feature in layer:
        geom = feature.GetGeometryRef()
        area = geom.GetArea()
        LayerArea = LayerArea + area
    
    return LayerArea


#--------------------------------------------------------------
def RandomInSituArea(shapefile, lc, code, crop, nbdraws, opath, workdir):
	"""
        This function creates 2 * nbdraws new shapefiles by selecting 50% of polygons of each crop class present for \n
        a learning file and the remaining 50% for a validation file
        ARGs:
		IN:
        - shapefile: the input shapefile
        - field: the name of the field in which selection will be based
        - nbdraws: the number of random selections wanted
        - opath: the output path
		OUTPUT:
        - 2 * nbdraws new shapefiles
        (From RandomSelectionInsitu_LV.py)
        """
	classes = []
	croptype = {}
        classesname = {}
	allFID = []
	nameshp = shapefile.split('.')
	namefile = '.'.join(nameshp[0:-1])
	namefile = namefile.split('/')
	
	driver = ogr.GetDriverByName("ESRI Shapefile")
	dataSource = driver.Open(shapefile, 0)
	layer = dataSource.GetLayer()
	
	# Select the crop classes and count the total features of cropland
	
	count = float(layer.GetFeatureCount())
	
	# Find the number of polygons by class
	for feature in layer:
		# Get the ID of the crop polygons and add them in a list
		pid = feature.GetFID()
		allFID.append(pid)
		# Get the list of codes of the crop class and add them in a list
		cl =  feature.GetField(field)
                ty =  feature.GetField(code)
                nm =  feature.GetField(lc)
		if cl not in classes:
			classes.append(cl)
                        croptype[cl] = ty
                        classesname[cl] = nm
   
        sortedclasses = sorted(classes)
        for cl in sortedclasses:
            print cl,croptype[cl],classesname[cl]

        file_Classes = open(workdir + "/Classes.txt","w") 
        file_Classes_CSV = open(workdir + "/Classes.csv","w") 
        file_Classes_CROP = open(workdir + "/Classes_CROP.txt","w") 
        file_Classes_CROP_CSV = open(workdir + "/Classes_CROP.csv","w") 
        file_Binary = open(workdir + "/Binary.txt","w") 
 
        for cl in sortedclasses:
            file_Classes.write("%s\n"%(cl)) 
            file_Classes_CSV.write("%s:%s\n"%(classesname[cl],cl)) 
            if croptype[cl] == 1:
                file_Classes_CROP.write("%s\n"%(cl)) 
                file_Classes_CROP_CSV.write("%s:%s\n"%(classesname[cl],cl)) 
            file_Binary.write("%s\n"%(croptype[cl])) 

        file_Classes_CSV.write("%s:%s\n"%("Nodecision",0))
        file_Classes_CROP_CSV.write("%s:%s\n"%("Nodecision",0))
 
        file_Classes.close() 
        file_Classes_CSV.close() 
        file_Classes_CROP.close() 
        file_Classes_CROP_CSV.close() 
        file_Binary.close() 
     
	#Crop classes are already selected, here after selects the individual classes of the crop classes
	for tirage in range(0,int(nbdraws)):
      		listallid = []
		for cl in classes:
			print "class " + str(cl)
			listid = []
            		listArea = []
			#Count the features of this class
			featureCount = float(layer.GetFeatureCount())
                	TotalArea = GetAreaOfPolygons(dataSource,field,cl)
            		#print cl
                	#print TotalArea
			#As we chose 50 percent, the total of features is divided by 2
			AreaThreshold = TotalArea *0.5  #-- (ADD Charlotte, for a percentage change : polbysel = int((featureCount*float(percentage))/100))			
			#Computes the proportion
			#Get all the IDs of the polygons of the class and add them in a list
			layer = dataSource.GetLayer()
			#Selects the polygons of the individual crop class
			layer.SetAttributeFilter(field+" = "+str(cl))
			for feat in layer:
				_id = feat.GetFID()
				listid.append(_id)
				geom = feat.GetGeometryRef()
				area = geom.GetArea()
                  		listArea.append(area)
				SEED = long(time.time()*256)
                                #SEED= 9091536107
				random.seed(SEED)
				random.shuffle(listid)
				random.seed(SEED)
				random.shuffle(listArea)
			Area = 0
			i = 0
			listToChoice = [] 
			while (Area < AreaThreshold):		
				listToChoice.append(listid[i])
				Area = Area + listArea[i]
				i = i + 1
			
			print "Total area " + str(TotalArea)
			print "Training area " + str(Area)
			print "Training polygons " + str(i-1)
	
			
			for fid in listToChoice:
				listallid.append(fid)
		listallid.sort()

		#Process to create the SQL query
		ch = ""
		listFid = []
		#Add the word FID=
		for fid in listallid:
			listFid.append("FID="+str(fid))
		#Add the word OR
		resultA = []
		for e in listFid:
			resultA.append(e)
			resultA.append(' OR ')
		resultA.pop()
		#Creates the SQL chain
		chA =  ''.join(resultA)
		#Select polygons with the SQL query
		layer.SetAttributeFilter(chA)
		print "namefile : ",namefile[-1]
		outShapefile = opath+"/Run_"+str(tirage)+".dir/"+namefile[-1]+"_seed"+str(tirage)+"_learn.shp"
		if os.path.isfile(outShapefile):
			continue
		#Create a new layer with the selection
		CreateNewLayer(layer, outShapefile)
		
		#Process to create the 2nd SQL request by choosing the polygons that were not choose before
		
		listValid = []
		for i in allFID:
			if i not in listallid:
				listValid.append(i)
        
		chV = ""
		#Add the word FID=
		listFidV = []
		for fid in listValid:
			listFidV.append("FID="+str(fid))
        
		resultV = []
		#Add the word OR
		for e in listFidV:
			resultV.append(e)
			resultV.append(' OR ')
		resultV.pop()
		
		chV =  ''.join(resultV)
		#Select polygons with the SQL query
		layer.SetAttributeFilter(chV)
		outShapefile2 = opath+"/Run_"+str(tirage)+".dir/"+namefile[-1]+"_seed"+str(tirage)+"_val.shp"
		CreateNewLayer(layer, lc, code, crop, outShapefile2)

#--------------------------------------------------------------
def CreateNewLayer(layer, lc, code, crop, outShapefile):
	"""
        This function creates a new shapefile
		ARGs:
        - layer: the input shapefile
        - outShapefile: the name of the output shapefile
		(From RandomSelectionInsitu_LV.py)
        """
	
	#Warning: used to S2AGRI data model, next line to change, modify name of attributs
	field_name_target = ['ID',lc, code, crop]
	#field_name_target = ['idParcelle','CodeS2agri','Nom_S2agri']
	outDriver = ogr.GetDriverByName("ESRI Shapefile")
	#if file already exists, delete it
	if os.path.exists(outShapefile):
		outDriver.DeleteDataSource(outShapefile)
	outDataSource = outDriver.CreateDataSource(outShapefile)
	out_lyr_name = os.path.splitext( os.path.split( outShapefile )[1] )[0]
	#Get the spatial reference of the input layer
	srsObj = layer.GetSpatialRef()
	#Creates the spatial reference of the output layer
	outLayer = outDataSource.CreateLayer( out_lyr_name, srsObj, geom_type=ogr.wkbMultiPolygon )
	# Add input Layer Fields to the output Layer if it is the one we want
	inLayerDefn = layer.GetLayerDefn()
	for i in range(0, inLayerDefn.GetFieldCount()):
		fieldDefn = inLayerDefn.GetFieldDefn(i)
		fieldName = fieldDefn.GetName()
		if fieldName not in field_name_target:
			continue
		outLayer.CreateField(fieldDefn)
	# Get the output Layer's Feature Definition
	outLayerDefn = outLayer.GetLayerDefn()
	cmpt = 1
	
	# Add features to the ouput Layer
	for inFeature in layer:
		# Create output Feature
		outFeature = ogr.Feature(outLayerDefn)
		# Add field values from input Layer
		for i in range(0, inLayerDefn.GetFieldCount()):
			fieldDefn = inLayerDefn.GetFieldDefn(i)
			fieldName = fieldDefn.GetName()
			if fieldName  in field_name_target:
				outFeature.SetField(inLayerDefn.GetFieldDefn(i).GetNameRef(),inFeature.GetField(i))
        
		# Set geometry as centroid
		geom = inFeature.GetGeometryRef()
		if geom == None:
			print cmpt
		outFeature.SetGeometry(geom.Clone())
		cmpt += 1
		# Add new feature to output Layer
		outLayer.CreateFeature(outFeature)	


if __name__ == "__main__":
 
    shpFile =  sys.argv[1]
    NbRun = int(sys.argv[2])
    OutputDir = sys.argv[3]
    WorkDir = sys.argv[4]
    try:
	lc = sys.argv[5]
	code = sys.argv[6]
	crop = sys.argv[7]
    except:
	lc = "LC"
	code = "CODE"
	crop = "CROP"
 

    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)

    for i in range(NbRun):
        if not os.path.exists(OutputDir + "/Run_%d.dir/"%(i) ):
            os.makedirs(OutputDir + "/Run_%d.dir/"%(i) ) 

   RandomInSituArea(shpFile, lc, code, crop, NbRun,OutputDir,WorkDir)



