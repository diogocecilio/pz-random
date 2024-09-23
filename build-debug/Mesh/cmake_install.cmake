# Install script for directory: /home/diogo/projects/pz-random/Mesh

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
    "/home/diogo/projects/pz-random/Mesh/TPZAgglomerateEl.h"
    "/home/diogo/projects/pz-random/Mesh/TPZCompElDisc.h"
    "/home/diogo/projects/pz-random/Mesh/TPZCompElLagrange.h"
    "/home/diogo/projects/pz-random/Mesh/TPZGeoElement.h"
    "/home/diogo/projects/pz-random/Mesh/TPZGeoElement.h.h"
    "/home/diogo/projects/pz-random/Mesh/TPZInterfaceEl.h"
    "/home/diogo/projects/pz-random/Mesh/TPZMultiphysicsInterfaceEl.h"
    "/home/diogo/projects/pz-random/Mesh/doxmesh.h"
    "/home/diogo/projects/pz-random/Mesh/pzcheckgeom.h"
    "/home/diogo/projects/pz-random/Mesh/pzcheckmesh.h"
    "/home/diogo/projects/pz-random/Mesh/pzcheckrestraint.h"
    "/home/diogo/projects/pz-random/Mesh/pzcmesh.h"
    "/home/diogo/projects/pz-random/Mesh/pzcompel.h"
    "/home/diogo/projects/pz-random/Mesh/pzcompelwithmem.h"
    "/home/diogo/projects/pz-random/Mesh/pzcondensedcompel.h"
    "/home/diogo/projects/pz-random/Mesh/pzconnect.h"
    "/home/diogo/projects/pz-random/Mesh/pzcreateapproxspace.h"
    "/home/diogo/projects/pz-random/Mesh/pzelchdiv.h"
    "/home/diogo/projects/pz-random/Mesh/pzelchdivbound2.h"
    "/home/diogo/projects/pz-random/Mesh/pzelctemp.h"
    "/home/diogo/projects/pz-random/Mesh/pzelementgroup.h"
    "/home/diogo/projects/pz-random/Mesh/pzelmat.h"
    "/home/diogo/projects/pz-random/Mesh/pzflowcmesh.h"
    "/home/diogo/projects/pz-random/Mesh/pzgeoel.h"
    "/home/diogo/projects/pz-random/Mesh/pzgeoelbc.h"
    "/home/diogo/projects/pz-random/Mesh/pzgeoelrefless.h"
    "/home/diogo/projects/pz-random/Mesh/pzgeoelrefless.h.h"
    "/home/diogo/projects/pz-random/Mesh/pzgeoelside.h"
    "/home/diogo/projects/pz-random/Mesh/pzgmesh.h"
    "/home/diogo/projects/pz-random/Mesh/pzgnode.h"
    "/home/diogo/projects/pz-random/Mesh/pzhdivfull.h"
    "/home/diogo/projects/pz-random/Mesh/pzhdivpressure.h"
    "/home/diogo/projects/pz-random/Mesh/pzhdivpressurebound.h"
    "/home/diogo/projects/pz-random/Mesh/pzintel.h"
    "/home/diogo/projects/pz-random/Mesh/pzinterpolationspace.h"
    "/home/diogo/projects/pz-random/Mesh/pzmeshid.h"
    "/home/diogo/projects/pz-random/Mesh/pzmultiphysicscompel.h"
    "/home/diogo/projects/pz-random/Mesh/pzmultiphysicselement.h"
    "/home/diogo/projects/pz-random/Mesh/pzreducedspace.h"
    "/home/diogo/projects/pz-random/Mesh/pzreferredcompel.h"
    "/home/diogo/projects/pz-random/Mesh/pzsubcmesh.h"
    "/home/diogo/projects/pz-random/Mesh/tpzagglomeratemesh.h"
    "/home/diogo/projects/pz-random/Mesh/tpzcompmeshreferred.h"
    "/home/diogo/projects/pz-random/Mesh/tpzgeoelmapped.h"
    "/home/diogo/projects/pz-random/Mesh/tpzgeoelrefpattern.h"
    "/home/diogo/projects/pz-random/Mesh/tpzgeoelrefpattern.h.h"
    )
endif()

