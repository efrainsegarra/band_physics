cmake_minimum_required (VERSION 2.6)
set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(single_photon_tdc single_photon_tdc.cpp )
target_link_libraries(single_photon_tdc ${ROOT_LIBRARIES} )

add_executable(single_photon single_photon.cpp )
target_link_libraries(single_photon ${ROOT_LIBRARIES} )

add_executable(corrected_tof corrected_tof.cpp )
target_link_libraries(corrected_tof ${ROOT_LIBRARIES} )
