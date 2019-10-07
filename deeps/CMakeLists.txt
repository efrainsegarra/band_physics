cmake_minimum_required (VERSION 2.6)
set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})


add_executable(w_plot w_plot.cpp )
target_link_libraries(w_plot ${ROOT_LIBRARIES} )
