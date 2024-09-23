# Install script for directory: /home/diogo/projects/pz-random/SpecialMaps

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
    "/home/diogo/projects/pz-random/SpecialMaps/TPZQuadSphere.h"
    "/home/diogo/projects/pz-random/SpecialMaps/TPZQuadTorus.h"
    "/home/diogo/projects/pz-random/SpecialMaps/TPZTriangleSphere.h"
    "/home/diogo/projects/pz-random/SpecialMaps/TPZTriangleTorus.h"
    "/home/diogo/projects/pz-random/SpecialMaps/TPZWavyLine.h"
    "/home/diogo/projects/pz-random/SpecialMaps/convtest.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzarc3d.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzblendnaca.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzchangeel.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzellipse3d.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzgeomid.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzmathtools.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzquadraticcube.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzquadraticline.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzquadraticprism.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzquadraticpyramid.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzquadraticquad.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzquadratictetra.h"
    "/home/diogo/projects/pz-random/SpecialMaps/tpzquadratictrig.h"
    )
endif()

