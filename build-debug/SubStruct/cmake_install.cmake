# Install script for directory: /home/diogo/projects/pz-random/SubStruct

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
    "/home/diogo/projects/pz-random/SubStruct/TPZTimeTemp.h"
    "/home/diogo/projects/pz-random/SubStruct/TPZfTime.h"
    "/home/diogo/projects/pz-random/SubStruct/doxsubstrmatrix.h"
    "/home/diogo/projects/pz-random/SubStruct/pzdohrstructmatrix.h"
    "/home/diogo/projects/pz-random/SubStruct/tpzdohrassemblelist.h"
    "/home/diogo/projects/pz-random/SubStruct/tpzdohrassembly.h"
    "/home/diogo/projects/pz-random/SubStruct/tpzdohrmatrix.h"
    "/home/diogo/projects/pz-random/SubStruct/tpzdohrprecond.h"
    "/home/diogo/projects/pz-random/SubStruct/tpzdohrsubstruct.h"
    "/home/diogo/projects/pz-random/SubStruct/tpzdohrsubstructCondense.h"
    "/home/diogo/projects/pz-random/SubStruct/tpzgensubstruct.h"
    "/home/diogo/projects/pz-random/SubStruct/tpzmatredstructmatrix.h"
    "/home/diogo/projects/pz-random/SubStruct/tpzpairstructmatrix.h"
    "/home/diogo/projects/pz-random/SubStruct/tpzparallelenviroment.h"
    )
endif()

