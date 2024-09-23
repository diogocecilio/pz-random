# Install script for directory: /home/diogo/projects/pz-random/Matrix

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
    "/home/diogo/projects/pz-random/Matrix/TPZCopySolve.h"
    "/home/diogo/projects/pz-random/Matrix/TPZPardisoSolver.h"
    "/home/diogo/projects/pz-random/Matrix/doxmatrix.h"
    "/home/diogo/projects/pz-random/Matrix/pzblock.h"
    "/home/diogo/projects/pz-random/Matrix/pzblockdiag.h"
    "/home/diogo/projects/pz-random/Matrix/pzbndmat.h"
    "/home/diogo/projects/pz-random/Matrix/pzdiffmatrix.h"
    "/home/diogo/projects/pz-random/Matrix/pzespmat.h"
    "/home/diogo/projects/pz-random/Matrix/pzfmatrix.h"
    "/home/diogo/projects/pz-random/Matrix/pzlink.h"
    "/home/diogo/projects/pz-random/Matrix/pzmatred.h"
    "/home/diogo/projects/pz-random/Matrix/pzmatrix.h"
    "/home/diogo/projects/pz-random/Matrix/pzmatrixid.h"
    "/home/diogo/projects/pz-random/Matrix/pzsbndmat.h"
    "/home/diogo/projects/pz-random/Matrix/pzseqsolver.h"
    "/home/diogo/projects/pz-random/Matrix/pzsespmat.h"
    "/home/diogo/projects/pz-random/Matrix/pzsfulmat.h"
    "/home/diogo/projects/pz-random/Matrix/pzshtmat.h"
    "/home/diogo/projects/pz-random/Matrix/pzskylmat.h"
    "/home/diogo/projects/pz-random/Matrix/pzskylmatpar.h"
    "/home/diogo/projects/pz-random/Matrix/pzskylnsymmat.h"
    "/home/diogo/projects/pz-random/Matrix/pzsolve.h"
    "/home/diogo/projects/pz-random/Matrix/pzspblockdiagpivot.h"
    "/home/diogo/projects/pz-random/Matrix/pzstencil.h"
    "/home/diogo/projects/pz-random/Matrix/pzstepsolver.h"
    "/home/diogo/projects/pz-random/Matrix/pzsysmp.h"
    "/home/diogo/projects/pz-random/Matrix/pzysmp.h"
    "/home/diogo/projects/pz-random/Matrix/tpzsparseblockdiagonal.h"
    "/home/diogo/projects/pz-random/Matrix/tpzverysparsematrix.h"
    )
endif()

