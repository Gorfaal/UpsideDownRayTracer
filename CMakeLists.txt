cmake_minimum_required(VERSION 2.8)

PROJECT(final)
SET(VTK_DIR /usr/local/lib/vtk-5.10)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(final final)

if(VTK_LIBRARIES)
  target_link_libraries(final ${VTK_LIBRARIES})
else()
  target_link_libraries(final vtkHybrid)
endif()

