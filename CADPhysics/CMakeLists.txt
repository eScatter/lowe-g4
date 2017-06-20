#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_CXX_FLAGS "-I$(SEM_BASE)/include")
add_definitions(-std=c++11)
set(BUILD_SHARED_LIBS ON)
project(CADPhysics)

find_package(Geant4 REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX)

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
                    3rdparty/csread
                    /home/annemarie/Downloads/hdf5-1.10.1/hdf5/include)
add_subdirectory(3rdparty/csread)
set(_hdf5_libs hdf5 hdf5_cpp)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/G4/src/*.cc ${PROJECT_SOURCE_DIR}/hdf5read/material.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/G4/include/*.hh ${PROJECT_SOURCE_DIR}/hdf5read/table/*.hh ${PROJECT_SOURCE_DIR}/hdf5read/units/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(CADPhysics ${sources} ${headers})
target_link_libraries(CADPhysics ${Geant4_LIBRARIES} )
target_link_libraries(CADPhysics csread /home/annemarie/Downloads/hdf5-1.10.1/hdf5/lib/libhdf5_cpp.so /home/annemarie/Downloads/hdf5-1.10.1/hdf5/lib/libhdf5.so -lz -ldl)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(DIRECTORY DESTINATION)
install(TARGETS CADPhysics DESTINATION lib)
install(DIRECTORY include DESTINATION .)
install(DIRECTORY mott DESTINATION share/CAD)
install(DIRECTORY df DESTINATION share/CAD)
