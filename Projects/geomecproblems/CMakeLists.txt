file(GLOB headers *.h)
file(GLOB sources *.cpp)

file (GLOB to_remove_from_sources main.cpp mainGUI.cpp poroelastoplastic.cpp mainPV.cpp)
list (REMOVE_ITEM sources ${to_remove_from_sources})

set(main1_SRCS main.cpp  ${sources} ${headers})

include_directories('src')
include_directories("/usr/include/eigen3")
#Plasticity works only with double type.
if ( STATE_TYPE STREQUAL "double" AND
     REAL_TYPE STREQUAL "double" AND 
     USING_FAD AND 
     BUILD_PLASTICITY_MATERIALS )


	  add_executable(geomecproblems ${main1_SRCS})
	  target_link_libraries(geomecproblems pz)



endif()

#include (CTestTestFile.cmake)
