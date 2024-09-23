# Install script for directory: /home/diogo/projects/pz-random

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/pzlib/include" TYPE FILE FILES "/home/diogo/projects/pz-random/build-debug/Common/config.h")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/diogo/projects/pz-random/build-debug/PerfUtil/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Util/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Common/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Save/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Integral/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/LinearSolvers/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Matrix/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Topology/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Geom/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/SpecialMaps/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Shape/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Refine/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/External/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Material/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Mesh/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Analysis/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Multigrid/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Post/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Frontal/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/StrMatrix/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/Pre/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/SubStruct/cmake_install.cmake")
  include("/home/diogo/projects/pz-random/build-debug/lib/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/diogo/projects/pz-random/build-debug/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
