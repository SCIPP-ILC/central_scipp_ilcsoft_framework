########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################

#Use a compiler that isn't outdated...
SET( CMAKE_CXX_COMPILER "/var/ilc/scipp/opt/gcc-5.3.0/build/bin/c++" )

# project name
PROJECT(scipp_ilcsoft_framework)


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

#Additional modules
set( CMAKE_SOURCE_PREFIX /var/ilc/scipp/central_scipp_ilcsoft_framework/modules_path )
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_PREFIX}")

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
FIND_PACKAGE( ROOT 5.34 REQUIRED)
FIND_PACKAGE( XercesC REQUIRED)
#FIND_PACKAGE( MarlinUtil 1.4 REQUIRED )
#FIND_PACKAGE( MarlinKinfit 0.0 REQUIRED )
#FIND_PACKAGE( GSL REQUIRED )

#FOREACH( pkg Marlin MarlinUtil MarlinKinfit GSL AIDA ROOT )
FOREACH( pkg Marlin ROOT XercesC)
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()


#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#   message(STATUS "dir='${dir}'")
#endforeach()


# optional package
#FIND_PACKAGE( AIDA )
#IF( AIDA_FOUND )
#    INCLUDE_DIRECTORIES( ${AIDA_INCLUDE_DIRS} )
#    LINK_LIBRARIES( ${AIDA_LIBRARIES} )
#    ADD_DEFINITIONS( "-DMARLIN_USE_AIDA" )
#    MESSAGE( STATUS "AIDA -- found" )
#ELSE()
#    MESSAGE( STATUS "AIDA -- not found" )
#ENDIF()



### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################

# definitions to pass to the compiler
# NOTE:I have altered these to use c++11 instead of ansi
ADD_DEFINITIONS( "-Wextra -std=c++11 -pedantic" ) #ADD_DEFINITIONS( "-Wall -ansi -pedantic" ) 
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ./src/core_processors/include ./src/processors/include ./src/base/util ./src/base/beamcal_recon )
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# source directories
AUX_SOURCE_DIRECTORY( ./src library_sources )
AUX_SOURCE_DIRECTORY( ./src/core_processors library_sources ) 
AUX_SOURCE_DIRECTORY( ./src/processors library_sources ) 
AUX_SOURCE_DIRECTORY( ./src/base/util library_sources )
AUX_SOURCE_DIRECTORY( ./src/base/beamcal_recon library_sources )

# add library
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib)

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

