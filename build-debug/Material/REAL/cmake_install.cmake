# Install script for directory: /home/diogo/projects/pz-random/Material/REAL

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
    "/home/diogo/projects/pz-random/Material/REAL/TPZConsLawTest.h"
    "/home/diogo/projects/pz-random/Material/REAL/TPZDarcyFlow.h"
    "/home/diogo/projects/pz-random/Material/REAL/TPZDarcyFlowIsoPerm.h"
    "/home/diogo/projects/pz-random/Material/REAL/TPZDiffusionConsLaw.h"
    "/home/diogo/projects/pz-random/Material/REAL/TPZEuler.h"
    "/home/diogo/projects/pz-random/Material/REAL/TPZFVHybrid.h"
    "/home/diogo/projects/pz-random/Material/REAL/TPZIsotropicPermeability.h"
    "/home/diogo/projects/pz-random/Material/REAL/TPZLinearConvection.h"
    "/home/diogo/projects/pz-random/Material/REAL/TPZMatElasticity2D.h"
    "/home/diogo/projects/pz-random/Material/REAL/TPZMulticamadaOrtho.h"
    "/home/diogo/projects/pz-random/Material/REAL/TPZPlacaOrthotropic.h"
    "/home/diogo/projects/pz-random/Material/REAL/eulerdif.h"
    "/home/diogo/projects/pz-random/Material/REAL/mixedpoisson.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzartdiff.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzausmflux.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzbctension.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzbiharmonic.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzblackoil2p3d.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzburger.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzconvectionproblem.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzcoupledtransportdarcy.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzcoupledtransportdarcyBC.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzdarcymatwithmem.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzdarcymem.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzelasAXImat.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzelasmat.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzelasmatwithmem.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzelast3d.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzelasthybrid.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzelasticmem.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzeuler.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzeulerconslaw.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzgradientflux.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzincnskeps.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzmaterialcoupling.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzmathyperelastic.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzmatmixedpoisson3d.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzmatorthotropic.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzmatplaca2.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzmattest.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzmattest3d.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzmatwithmem.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzmultiphase.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzmultplaca.h"
    "/home/diogo/projects/pz-random/Material/REAL/pznlmat1d.h"
    "/home/diogo/projects/pz-random/Material/REAL/pznlmat1drotatedengstrain.h"
    "/home/diogo/projects/pz-random/Material/REAL/pznonlinbiharmonic.h"
    "/home/diogo/projects/pz-random/Material/REAL/pznonlinearpoisson3d.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzplaca.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzpoisson3dreferred.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzspacetimerichardseq.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzthermicelast3d.h"
    "/home/diogo/projects/pz-random/Material/REAL/pzviscoelastic.h"
    "/home/diogo/projects/pz-random/Material/REAL/swelling.h"
    "/home/diogo/projects/pz-random/Material/REAL/tpzmultcamada.h"
    )
endif()

