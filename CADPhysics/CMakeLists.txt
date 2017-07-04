#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_CXX_FLAGS "-I$(SEM_BASE)/include")
add_definitions(-std=c++11)
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
                    ${PROJECT_SOURCE_DIR}/G4/include
                    3rdparty/csread)
add_subdirectory(3rdparty/csread)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/G4/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/G4/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(CADPhysics ${sources} ${headers})
target_link_libraries(CADPhysics ${Geant4_LIBRARIES})
set_property(TARGET csread PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(CADPhysics csread)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(DIRECTORY DESTINATION)
install(TARGETS CADPhysics DESTINATION lib)
install(DIRECTORY include DESTINATION .)
