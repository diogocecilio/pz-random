# Install script for directory: /home/diogo/projects/pz-random/Pre

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
    "/home/diogo/projects/pz-random/Pre/TPZExtendGridDimension.h"
    "/home/diogo/projects/pz-random/Pre/TPZGMSHReadMesh.h"
    "/home/diogo/projects/pz-random/Pre/TPZGenSpecialGrid.h"
    "/home/diogo/projects/pz-random/Pre/TPZMHMeshControl.h"
    "/home/diogo/projects/pz-random/Pre/TPZReadGIDGrid.h"
    "/home/diogo/projects/pz-random/Pre/doxpre.h"
    "/home/diogo/projects/pz-random/Pre/pzbuildmultiphysicsmesh.h"
    "/home/diogo/projects/pz-random/Pre/pzdatafi.h"
    "/home/diogo/projects/pz-random/Pre/pzgengrid.h"
    "/home/diogo/projects/pz-random/Pre/pzhyperplane.h"
    "/home/diogo/projects/pz-random/Pre/pzidentifyrefpattern.h"
    "/home/diogo/projects/pz-random/Pre/pzpargrid.h"
    "/home/diogo/projects/pz-random/Pre/pzreadmesh.h"
    "/home/diogo/projects/pz-random/Pre/pzreadmeshhr.h"
    "/home/diogo/projects/pz-random/Pre/pzreadtetgen.h"
    "/home/diogo/projects/pz-random/Pre/tpzhierarquicalgrid.h"
    )
endif()

