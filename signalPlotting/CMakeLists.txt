
cmake_minimum_required (VERSION 2.6)
set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)
include_directories(${INCLUDES})


add_executable(2d_reach 2d_reach.cpp )
target_link_libraries(2d_reach ${ROOT_LIBRARIES} )
