# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( GSoC_ )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  create_single_source_cgal_program( "benchmark_discrete_distrib.cpp" )
  create_single_source_cgal_program( "mesh-3D-sampling-test.cpp" )
  create_single_source_cgal_program( "mesh-3D-sampling-test-improved.cpp" )
  create_single_source_cgal_program( "points_in_3D_mesh.cpp" )
  create_single_source_cgal_program( "points_in_3_IFs.cpp" )
  create_single_source_cgal_program( "points_in_triangle_2.cpp" )
  create_single_source_cgal_program( "points_in_triangle_3.cpp" )
  create_single_source_cgal_program( "points_in_triangle.cpp" )
  create_single_source_cgal_program( "random_disc_2.cpp" )
  create_single_source_cgal_program( "random_grid.cpp" )
  create_single_source_cgal_program( "surface_mesh_sampling_test.cpp" )
  create_single_source_cgal_program( "triangle_3_test.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

