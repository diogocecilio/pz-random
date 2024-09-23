# Install script for directory: /home/diogo/projects/pz-random/StrMatrix

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
    "/home/diogo/projects/pz-random/StrMatrix/KLStrMatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/TPBSpStructMatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/TPZFrontStructMatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/TPZParFrontStructMatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/TPZParSkylineStructMatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/TPZSkylineNSymStructMatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/TPZSpStructMatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/doxstrmatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/pzbdstrmatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/pzbstrmatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/pzequationfilter.h"
    "/home/diogo/projects/pz-random/StrMatrix/pzfstrmatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/pzsbstrmatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/pzskylstrmatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/pzstrmatrix.h"
    "/home/diogo/projects/pz-random/StrMatrix/pzstrmatrixcs.h"
    "/home/diogo/projects/pz-random/StrMatrix/pzstrmatrixgc.h"
    "/home/diogo/projects/pz-random/StrMatrix/pzstrmatrixot.h"
    "/home/diogo/projects/pz-random/StrMatrix/pzstrmatrixst.h"
    "/home/diogo/projects/pz-random/StrMatrix/pzstrmatrixtbb.h"
    "/home/diogo/projects/pz-random/StrMatrix/tpzsparseblockdiagonalstructmatrix.h"
    )
endif()

