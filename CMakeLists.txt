# CMakeLists for Deuteron DIS project

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

project (analysis)

add_subdirectory(backgroundSubtraction)
add_subdirectory(eventMixing)
add_subdirectory(acceptance)
add_subdirectory(signalPlotting)
add_subdirectory(photonPeak)
add_subdirectory(crossSection)
add_subdirectory(deeps)

