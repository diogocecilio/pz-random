# Install script for directory: /home/diogo/projects/pz-random/Material

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
    "/home/diogo/projects/pz-random/Material/TPZElast3Dnlinear.h"
    "/home/diogo/projects/pz-random/Material/TPZLagrangeMultiplier.h"
    "/home/diogo/projects/pz-random/Material/TPZLinearConvecDiff.h"
    "/home/diogo/projects/pz-random/Material/TPZMatDualHybridPoisson.h"
    "/home/diogo/projects/pz-random/Material/TPZMatLaplacian.h"
    "/home/diogo/projects/pz-random/Material/TPZReynoldsFlow.h"
    "/home/diogo/projects/pz-random/Material/doxmaterial.h"
    "/home/diogo/projects/pz-random/Material/pzbndcond.h"
    "/home/diogo/projects/pz-random/Material/pzconslaw.h"
    "/home/diogo/projects/pz-random/Material/pzdiscgal.h"
    "/home/diogo/projects/pz-random/Material/pzl2projection.h"
    "/home/diogo/projects/pz-random/Material/pzmat1dlin.h"
    "/home/diogo/projects/pz-random/Material/pzmat2dlin.h"
    "/home/diogo/projects/pz-random/Material/pzmaterial.h"
    "/home/diogo/projects/pz-random/Material/pzmaterialdata.h"
    "/home/diogo/projects/pz-random/Material/pzmaterialid.h"
    "/home/diogo/projects/pz-random/Material/pzpoisson3d.h"
    "/home/diogo/projects/pz-random/Material/pztransientmat.h"
    "/home/diogo/projects/pz-random/Material/pzuncoupledmultiphysics.h"
    "/home/diogo/projects/pz-random/Material/tpzoutofrange.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/diogo/projects/pz-random/build-debug/Material/REAL/cmake_install.cmake")

endif()

