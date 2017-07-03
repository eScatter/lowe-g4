#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_CXX_FLAGS "-I/usr/local/include -I/usr/local/include/G4 -I$(SEM_BASE)/include")
project(SEM_4.10.00)

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
                    ${PROJECT_SOURCE_DIR}/include
                    /home/annemarie/Downloads/hdf5-1.10.1/hdf5/include)
set(_hdf5_libs hdf5 hdf5_cpp)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(SEM_4.10.00 SEM_4.10.00.cc ${sources})
target_link_libraries(SEM_4.10.00 /usr/local/lib/libCADPhysics.so ${Geant4_LIBRARIES}
                      /usr/local/lib/libcsread.so
                      /home/annemarie/Downloads/hdf5-1.10.1/hdf5/lib/libhdf5_cpp.so
                      /home/annemarie/Downloads/hdf5-1.10.1/hdf5/lib/libhdf5.so -lz -ldl)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS SEM_4.10.00 DESTINATION bin)
