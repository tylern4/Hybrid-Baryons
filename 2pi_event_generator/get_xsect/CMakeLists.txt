cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}" )
include_directories(../)
include_directories(../interpol)

add_library(getXsect STATIC
  get_xsect_14_18_lowq2_fit.cxx
  get_xsect_25_30.cxx
  get_xsect_fedotov.cxx
  get_xsect_gol2.cxx
  get_xsect_golovach.cxx
  get_xsect_near_threshold.cxx
  get_xsect_q2_130_w_gt_18_lt_21.cxx
  get_xsect_q2_13_wgt_3.cxx
  get_xsect_rip_fed_join.cxx
  get_xsect_ripani.cxx
  get_xsect_scale_gol_18_25.cxx
  )
