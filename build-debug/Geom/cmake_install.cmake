# Install script for directory: /home/diogo/projects/pz-random/Geom

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/pzlib/include" TYPE FILE FILES
    "/home/diogo/projects/pz-random/Geom/TPZGeoCube.h"
    "/home/diogo/projects/pz-random/Geom/TPZGeoExtend.h"
    "/home/diogo/projects/pz-random/Geom/TPZGeoLinear.h"
    "/home/diogo/projects/pz-random/Geom/doxgeometry.h"
    "/home/diogo/projects/pz-random/Geom/pzgeopoint.h"
    "/home/diogo/projects/pz-random/Geom/pzgeoprism.h"
    "/home/diogo/projects/pz-random/Geom/pzgeopyramid.h"
    "/home/diogo/projects/pz-random/Geom/pzgeoquad.h"
    "/home/diogo/projects/pz-random/Geom/pzgeotetrahedra.h"
    "/home/diogo/projects/pz-random/Geom/pzgeotriangle.h"
    "/home/diogo/projects/pz-random/Geom/pznoderep.h"
    "/home/diogo/projects/pz-random/Geom/pznoderep.h.h"
    "/home/diogo/projects/pz-random/Geom/pzshapeextend.h"
    "/home/diogo/projects/pz-random/Geom/tpzgeoblend.h"
    )
endif()

