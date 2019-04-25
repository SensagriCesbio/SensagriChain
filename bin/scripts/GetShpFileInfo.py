#!/usr/bin/python
import os, sys
from collections import Counter


from osgeo import ogr, osr
try:
    from osgeo import gdal
except ImportError:
    import gdal


#--------------------------------------------------------------------
#----This function returns the list of elements of the field "CODE" of the shape file
#--------------------------------------------------------------------

def getFieldElement(shape,driverName="ESRI Shapefile",field = "CODE",mode = "all",elemType = "int"):
    """
        IN :
        shape [string] : shape to compute
        driverName [string] : ogr driver to read the shape
        field [string] : data's field
        mode [string] : "all" or "unique"
        OUT :
        [list] containing all/unique element in shape's field
        
        Example :
        getFieldElement("./MyShape.sqlite","SQLite","CODE",mode = "all")
        >> [1,2,2,2,2,3,4]
        getFieldElement("./MyShape.sqlite","SQLite","CODE",mode = "unique")
        >> [1,2,3,4]
        """
    def getElem(elem,elemType):
        if elemType == "int" : return int(elem)
        elif elemType == "str" : return str(elem)
        else:
            raise Exception("elemType must be 'int' or 'str'")
    driver = ogr.GetDriverByName(driverName)
    dataSource = driver.Open(shape, 0)
    layer = dataSource.GetLayer()
    if mode == "all" : return [ getElem(currentFeat.GetField(field),elemType) for currentFeat in layer]
    elif mode == "unique" : return list(set([ getElem(currentFeat.GetField(field),elemType) for currentFeat in layer]))
    else:
        raise Exception("mode parameter must be 'all' or 'unique'")



#--------------------------------------------------------------------
#---This function returns the land cover dictionary of the shapeFile
#-- It allows us to know that the code 1313 corresponds to class Maize
#--------------------------------------------------------------------

def getLandCoverDictionary(shapeFile):

    LegendCodeList = getFieldElement(shapeFile,"ESRI Shapefile","CODE",mode = "unique")
    LCdictionary = dict([(key, []) for key in LegendCodeList])
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapeFile , 0)

    layer = dataSource.GetLayer()

    for feature in layer:
        idCode = feature.GetField("CODE")
        data = LCdictionary.get(idCode, "")

        if data == []:
            LCdictionary[idCode] = feature.GetField("LC")

    return LCdictionary

#--------------------------------------------------------------------
#---This function returns a binary dictionary of the shapeFile
#-- It allows us to know if the code 1313 corresponds to a crop class
#--------------------------------------------------------------------

def getCropCoverDictionary(shapeFile):
    
    LegendCodeList = getFieldElement(shapeFile,"ESRI Shapefile","CODE",mode = "unique")
    Dictionary = dict([(key, []) for key in LegendCodeList])
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapeFile , 0)
    
    layer = dataSource.GetLayer()
    
    for feature in layer:
        idCode = feature.GetField("CODE")
        data = Dictionary.get(idCode, "")
        
        if data == []:
            Dictionary[idCode] = feature.GetField("CROP")

    return Dictionary

def getClassNumber(shapeFile):
    
    LegendCodeList = getFieldElement(shapeFile,"ESRI Shapefile","CODE",mode = "unique")
    Dictionary = dict([(key, []) for key in LegendCodeList])
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapeFile , 0)
    
    layer = dataSource.GetLayer()
   
    classList = [  feature.GetField("CODE") for feature in layer]
    
    LabelCounter = Counter(classList)
    print LabelCounter

#--------------------------------------------------------------------

if __name__ == "__main__":
    # Get the input Layer
    inShapefile = sys.argv[1]
    
    CropDictionary =  getCropCoverDictionary(inShapefile)
    CoverDictionary = getLandCoverDictionary(inShapefile)

    print CropDictionary
    print CoverDictionary      

    getClassNumber(inShapefile) 
