#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(Hadr06)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
# NEW ADD 2021 10 14
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/channeling/include)
# End of new add

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
#include_directories(${PROJECT_SOURCE_DIR}/include 
#                    ${Geant4_INCLUDE_DIR})
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# NEW ADD 2021 10 14
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc 
                  ${PROJECT_SOURCE_DIR}/channeling/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh 
                  ${PROJECT_SOURCE_DIR}/channeling/include/*.hh)
                  
# Install data 
#
file(GLOB _data RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/data/*.txt)
foreach(_dat ${_data})
     configure_file(
	${PROJECT_SOURCE_DIR}/${_dat}
     	${PROJECT_BINARY_DIR}/${_dat}
	COPYONLY
	)
endforeach()

file(GLOB _macs RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/mac/*.mac)
foreach(_mac ${_macs})
     configure_file(
	${PROJECT_SOURCE_DIR}/${_mac}
     	${PROJECT_BINARY_DIR}/${_mac}
	COPYONLY
	)
endforeach()                  
                  
                  
                  
# End of new add
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(Hadr06 Hadr06.cc ${sources} ${headers})
target_link_libraries(Hadr06 ${Geant4_LIBRARIES} )
# NEW ADD 2021 10 14
message(${Geant4_LIBRARIES})
# target_link_libraries(channeling ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Hadr06. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(Hadr04_SCRIPTS
    vis.mac
  )

foreach(_script ${Hadr04_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS Hadr06 DESTINATION bin)

