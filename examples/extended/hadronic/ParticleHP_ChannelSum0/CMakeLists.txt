# Standalone repro for G4ParticleHPChannel isotope sum==0 (PR discussion).
# Not registered in the parent hadronic examples CMakeLists — build against an install:
#
#   cmake -S . -B build -DGeant4_DIR=$G4INSTALL/lib/cmake/Geant4
#   cmake --build build
#   ./build/ChannelSum0

cmake_minimum_required(VERSION 3.16...3.27)
project(ParticleHP_ChannelSum0)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(ChannelSum0 ChannelSum0.cc)
target_link_libraries(ChannelSum0 ${Geant4_LIBRARIES})

add_executable(ChannelSum0_scan ChannelSum0_scan.cc)
target_link_libraries(ChannelSum0_scan ${Geant4_LIBRARIES})
