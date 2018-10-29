cmake_minimum_required(VERSION 3.5)
project(2piEG)
set(VERSION 1.0)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net EG Hist)
include(${ROOT_USE_FILE})

set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${PROJECT_SOURCE_DIR}/cmake)

add_subdirectory(2pi_event_generator)
