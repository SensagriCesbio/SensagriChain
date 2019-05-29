cmake_minimum_required(VERSION 2.8)  
 
SET(CMAKE_CXX_FLAGS -std=c++11)
 
FIND_PACKAGE(OTB PATHS /usr/p18004/larnaud/otb-6.2.0-linux64/lib/cmake/otb-6.2/)  

IF(OTB_FOUND)  
  INCLUDE(${OTB_USE_FILE})  
ELSE(OTB_FOUND)  
  MESSAGE(FATAL_ERROR  
      "Cannot build OTB project without OTB.  Please set OTB_DIR.")  
ENDIF(OTB_FOUND)

option(OTB_USE_ITK "Enable module OpenCV in OTB" ON)
option(OTB_USE_OPENCV "Enable module OpenCV in OTB" ON)

SET(OTB_DIR /usr/p18004/larnaud/otb-6.2.0-linux64/lib/cmake/OTB-6.2)
SET(ITK_DIR /usr/p18004/larnaud/otb-6.2.0-linux64/lib/cmake/ITK-4.12)

otb_create_application(NAME  RFClassification
  SOURCES        otbRFClassification.cxx CvRTreesWrapper.cxx
  LINK_LIBRARIES ${OTB_LIBRARIES}) 
