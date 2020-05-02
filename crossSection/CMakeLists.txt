cmake_minimum_required (VERSION 2.6)
set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(count_events count_events.cpp )
target_link_libraries(count_events ${ROOT_LIBRARIES} )
