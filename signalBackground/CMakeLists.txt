cmake_minimum_required (VERSION 2.6)
set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(signal_background signal_background.cpp )
target_link_libraries(signal_background ${ROOT_LIBRARIES} )
