# Install script for directory: /home/diogo/projects/pz-random/Post

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
    "/home/diogo/projects/pz-random/Post/TPZProjectEllipse.h"
    "/home/diogo/projects/pz-random/Post/TPZVTKGeoMesh.h"
    "/home/diogo/projects/pz-random/Post/doxpost.h"
    "/home/diogo/projects/pz-random/Post/pzcompelpostproc.h"
    "/home/diogo/projects/pz-random/Post/pzdxmesh.h"
    "/home/diogo/projects/pz-random/Post/pzgradientreconstruction.h"
    "/home/diogo/projects/pz-random/Post/pzgraphel.h"
    "/home/diogo/projects/pz-random/Post/pzgraphel1d.h"
    "/home/diogo/projects/pz-random/Post/pzgraphel1dd.h"
    "/home/diogo/projects/pz-random/Post/pzgraphelq2d.h"
    "/home/diogo/projects/pz-random/Post/pzgraphelq2dd.h"
    "/home/diogo/projects/pz-random/Post/pzgraphelq3dd.h"
    "/home/diogo/projects/pz-random/Post/pzgraphmesh.h"
    "/home/diogo/projects/pz-random/Post/pzgraphnode.h"
    "/home/diogo/projects/pz-random/Post/pzmvmesh.h"
    "/home/diogo/projects/pz-random/Post/pzpostprocanalysis.h"
    "/home/diogo/projects/pz-random/Post/pzpostprocmat.h"
    "/home/diogo/projects/pz-random/Post/pztrigraph.h"
    "/home/diogo/projects/pz-random/Post/pztrigraphd.h"
    "/home/diogo/projects/pz-random/Post/pzv3dmesh.h"
    "/home/diogo/projects/pz-random/Post/pzvisualmatrix.h"
    "/home/diogo/projects/pz-random/Post/pzvtkmesh.h"
    "/home/diogo/projects/pz-random/Post/tpzgraphelprismmapped.h"
    "/home/diogo/projects/pz-random/Post/tpzgraphelpyramidmapped.h"
    "/home/diogo/projects/pz-random/Post/tpzgraphelt2dmapped.h"
    "/home/diogo/projects/pz-random/Post/tpzgraphelt3d.h"
    )
endif()

