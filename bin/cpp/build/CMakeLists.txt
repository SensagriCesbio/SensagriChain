CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(SenSAgri)

SET(CMKE_CXX_FLAGS -std=c++14)

FIND_PACKAGE(OTB PATHS /users/p18004/larnaud/OTB-6.2.0-Linux64/lib/cmake/OTB-6.2/)

IF(OTB_FOUND)
    INCLUDE(${OTB_USE_FILE})
ELSE(OTB_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build OTB project without OTB")
ENDIF(OTB_FOUND)

OPTION(OTB_USE_ITK "Enable module OpenCV in OTB" ON)
OPTION(OTB_USE_OPENCV "Enable module OpenCV in OTB" ON)

SET(OTB_DIR /users/p18004/larnaud/OTB-6.2.0-Linux64/lib/cmake/OTB-6.2)
SET(ITK_DIR /users/p18004/larnaud/OTB-6.2.0-Linux64/lib/cmake/ITK-4.12)

#INCLUDE_DIRECTORIES(${OTB_USE_FILE} "${CMAKE_CURRENT_SOURCE_DIR}/includes" "${CMAKE_CURRENT_SOURCE_DIR}/src"   )
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/includes" "${CMAKE_CURRENT_SOURCE_DIR}/src")

FILE(GLOB_RECURSE include_H  ${CMAKE_CURRENT_SOURCE_DIR}/includes/*.hpp )
INSTALL(FILES ${include_H} DESTINATION includes)

FILE(GLOB_RECURSE project_sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx )

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Executables/)

file( GLOB mainPrograms RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/MainPrograms" "${CMAKE_CURRENT_SOURCE_DIR}/MainPrograms/*.cxx" )

foreach( testsourcefile ${mainPrograms} )
    # I used a simple string replace, to cut off .cpp.
    string( REPLACE ".cxx" "" testname ${testsourcefile} )
    add_executable(${testname} ${CMAKE_CURRENT_SOURCE_DIR}/MainPrograms/${testsourcefile}  ${project_sources} )
    # Make sure YourLib is linked to each app
    target_link_libraries( ${testname} ${OTB_LIBRARIES} )
endforeach( testsourcefile ${mainPrograms} )
