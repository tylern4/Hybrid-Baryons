cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}" )
include_directories(../get_xsect)
include_directories(../)

add_library(interpol STATIC
  interpol.cxx
  interpol_fedotov.cxx
  interpol_fedotov_thresh.cxx
  interpol_gol2.cxx
  interpol_golovach.cxx
  interpol_int.cxx
  interpol_phot_wgt_3.cxx
  interpol_q2_13_wgt_3.cxx
  interpol_rip2.cxx
  interpol_rip3.cxx
  )
