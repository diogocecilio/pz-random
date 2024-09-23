# Install script for directory: /home/diogo/projects/pz-random/Util

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
    "/home/diogo/projects/pz-random/Util/TPZSemaphore.h"
    "/home/diogo/projects/pz-random/Util/TPZThreadTools.h"
    "/home/diogo/projects/pz-random/Util/TPZTimer.h"
    "/home/diogo/projects/pz-random/Util/checkconv.h"
    "/home/diogo/projects/pz-random/Util/doxutil.h"
    "/home/diogo/projects/pz-random/Util/pzadmchunk.h"
    "/home/diogo/projects/pz-random/Util/pzadmchunkthreadsafe.h"
    "/home/diogo/projects/pz-random/Util/pzaxestools.h"
    "/home/diogo/projects/pz-random/Util/pzcheckconsistency.h"
    "/home/diogo/projects/pz-random/Util/pzchunk.h"
    "/home/diogo/projects/pz-random/Util/pzfunction.h"
    "/home/diogo/projects/pz-random/Util/pzgradient.h"
    "/home/diogo/projects/pz-random/Util/pzline.h"
    "/home/diogo/projects/pz-random/Util/pzlog.h"
    "/home/diogo/projects/pz-random/Util/pzmanvector.h"
    "/home/diogo/projects/pz-random/Util/pznuma.h"
    "/home/diogo/projects/pz-random/Util/pznumeric.h"
    "/home/diogo/projects/pz-random/Util/pzpix.h"
    "/home/diogo/projects/pz-random/Util/pzplane.h"
    "/home/diogo/projects/pz-random/Util/pzpolynomial.h"
    "/home/diogo/projects/pz-random/Util/pzstack.h"
    "/home/diogo/projects/pz-random/Util/pzstring.h"
    "/home/diogo/projects/pz-random/Util/pzvec.h"
    "/home/diogo/projects/pz-random/Util/pzvec_extras.h"
    "/home/diogo/projects/pz-random/Util/tpzautopointer.h"
    "/home/diogo/projects/pz-random/Util/tpzpagemigrationmanager.h"
    "/home/diogo/projects/pz-random/Util/tpzpermutation.h"
    )
endif()

