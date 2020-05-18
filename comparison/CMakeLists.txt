
cmake_minimum_required (VERSION 2.6)
set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)
include_directories(${INCLUDES})


add_executable(data_sim_v0 data_sim_v0.cpp )
target_link_libraries(data_sim_v0 ${ROOT_LIBRARIES} )
