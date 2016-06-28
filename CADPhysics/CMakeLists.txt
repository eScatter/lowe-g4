#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_CXX_FLAGS "-I$(SEM_BASE)/include")
set(BUILD_SHARED_LIBS ON)
project(CADPhysics)

find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${PROJECT_SOURCE_DIR}/G4/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/G4/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/G4/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(CADPhysics ${sources} ${headers})
target_link_libraries(CADPhysics ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(DIRECTORY DESTINATION)
install(TARGETS CADPhysics DESTINATION lib)
install(DIRECTORY include DESTINATION .)
install(DIRECTORY mott DESTINATION share/CAD)
install(DIRECTORY df DESTINATION share/CAD)

