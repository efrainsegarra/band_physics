
cmake_minimum_required (VERSION 2.6)
set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(testSubtraction testSubtraction.cpp )
target_link_libraries(testSubtraction ${ROOT_LIBRARIES} )
