cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}" )
include_directories(get_xsect)
include_directories(interpol)

add_library(2pilib STATIC
  anti_rot.cxx
  fermi_anti_rot.cxx
  fermi_bonn.cxx
  fermi_rot.cxx
  global.cxx
  hist_def.cxx
  hist_fill.cxx
  hist_write.cxx
  inp_file_read.cxx
  out_file_close.cxx
  out_file_fill.cxx
  out_file_open.cxx
  radcorr.cxx
  read_fit_param_files.cxx
  read_xsect_files.cxx
  rot.cxx
  )
