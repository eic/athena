cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# CMP0074: find_package() uses <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

PROJECT(athena
  VERSION 0.0.1
  LANGUAGES CXX
  DESCRIPTION "A template dd4hep+acts detector"
  )

option(IRT_AUXFILE "Build auxiliary config file for dRICH IRT" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")

find_package(DD4hep REQUIRED COMPONENTS DDCore DDG4)
find_package(ActsDD4hep)
if(ActsDD4hep_FOUND)
  add_compile_definitions(USE_ACTSDD4HEP)
  set(ActsDD4hep ActsDD4hep::ActsDD4hep)
else()
  find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep)
  set(ActsDD4hep ActsCore ActsPluginDD4hep)
endif()
find_package(fmt REQUIRED)

if(IRT_AUXFILE)
  add_compile_definitions(IRT_AUXFILE)
  find_package(IRT REQUIRED)
endif()

#-----------------------------------------------------------------------------------
set(a_lib_name athena)

dd4hep_configure_output()

dd4hep_add_plugin(${a_lib_name} SOURCES 
  src/BarrelBarDetectorWithSideFrame_geo.cpp
  src/BarrelCalorimeter_geo.cpp
  src/BarrelCalorimeterInterlayers_geo.cpp
  src/BarrelTrackerWithFrame_geo.cpp
  src/CompositeTracker_geo.cpp
  src/CylinderTrackerBarrel_geo.cpp
  src/DIRC_geo.cpp
  src/DRICH_geo.cpp
  src/PFRICH_geo.cpp
  src/FileLoader.cpp
  src/FieldMapBrBz.cpp
  src/GaseousRICH_geo.cpp
  src/GeometryHelpers.cpp
  src/HomogeneousCalorimeter_geo.cpp
  src/HybridCalorimeter_geo.cpp
  src/MRich_geo.cpp
  src/PolyhedraEndcapCalorimeter2_geo.cpp
  src/ScFiCalorimeter_geo.cpp
  src/ShashlikCalorimeter_geo.cpp
  src/SimpleDiskDetector_geo.cpp
  src/SolenoidCoil_geo.cpp
  src/TrapEndcapTracker_geo.cpp
  USES ${ActsDD4hep}
  )
target_link_libraries(${a_lib_name}
  PUBLIC DD4hep::DDCore  DD4hep::DDRec fmt::fmt
  )
if(IRT_AUXFILE)
  target_link_libraries(${a_lib_name} PUBLIC IRT)
endif()

#-----------------------------------------------------------------------------------
# Install the detector description files.
install(DIRECTORY compact/
    DESTINATION share/${PROJECT_NAME}/${a_lib_name}
    FILES_MATCHING PATTERN "*.xml"
    )
install(DIRECTORY compact/
    DESTINATION share/${PROJECT_NAME}/compact
    FILES_MATCHING PATTERN "*.xml"
    )
install(FILES athena.xml
    DESTINATION share/${PROJECT_NAME}/
    )
#-----------------------------------------------------------------------------------
# Install the detector calibration files.
install(DIRECTORY calibrations/
    DESTINATION share/${PROJECT_NAME}/calibrations
    )
