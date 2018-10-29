cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}" )
include_directories(.)
include_directories(get_xsect)
include_directories(interpol)

add_executable(twopeg.exe main_prog.cxx
  global.cxx
  hist_def.cxx
  inp_file_read.cxx
  read_xsect_files.cxx
  read_fit_param_files.cxx
  out_file_open.cxx
  out_file_fill.cxx
  out_file_close.cxx
  get_xsect/get_xsect_ripani.cxx
  get_xsect/get_xsect_fedotov.cxx
  get_xsect/get_xsect_q2_130_w_gt_18_lt_21.cxx
  get_xsect/get_xsect_rip_fed_join.cxx
  get_xsect/get_xsect_near_threshold.cxx
  get_xsect/get_xsect_golovach.cxx
  get_xsect/get_xsect_14_18_lowq2_fit.cxx
  get_xsect/get_xsect_scale_gol_18_25.cxx
  anti_rot.cxx
  rot.cxx
  interpol/interpol.cxx
  interpol/interpol_rip2.cxx
  interpol/interpol_golovach.cxx
  interpol/interpol_fedotov.cxx
  interpol/interpol_fedotov_thresh.cxx
  interpol/interpol_rip3.cxx
  interpol/interpol_q2_13_wgt_3.cxx
  interpol/interpol_phot_wgt_3.cxx
  interpol/interpol_gol2.cxx
  get_xsect/get_xsect_gol2.cxx
  get_xsect/get_xsect_25_30.cxx
  interpol/interpol_int.cxx
  radcorr.cxx
  fermi_bonn.cxx
  fermi_rot.cxx
  fermi_anti_rot.cxx
  get_xsect/get_xsect_q2_13_wgt_3.cxx
  )

  target_link_libraries(twopeg.exe ${ROOT_LIBRARIES})
