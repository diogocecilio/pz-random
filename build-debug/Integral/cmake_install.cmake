# Install script for directory: /home/diogo/projects/pz-random/Integral

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
    "/home/diogo/projects/pz-random/Integral/adapt.h"
    "/home/diogo/projects/pz-random/Integral/doxintegral.h"
    "/home/diogo/projects/pz-random/Integral/pzquad.h"
    "/home/diogo/projects/pz-random/Integral/tpzgaussrule.h"
    "/home/diogo/projects/pz-random/Integral/tpzintpoints.h"
    "/home/diogo/projects/pz-random/Integral/tpzintrulelist.h"
    "/home/diogo/projects/pz-random/Integral/tpzintrulep3d.h"
    "/home/diogo/projects/pz-random/Integral/tpzintrulet.h"
    "/home/diogo/projects/pz-random/Integral/tpzintrulet3d.h"
    "/home/diogo/projects/pz-random/Integral/tpzprinteg.h"
    )
endif()

