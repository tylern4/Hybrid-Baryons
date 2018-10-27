cmake_minimum_required(VERSION 3.5)
project(2piEG)
set(VERSION 1.0)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
message(STATUS "ROOTSYS: " $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net EG)
include(${ROOT_USE_FILE})

set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${PROJECT_SOURCE_DIR}/cmake)


add_subdirectory(2pi_event_generator/get_xsect)
add_subdirectory(2pi_event_generator/interpol)
add_subdirectory(2pi_event_generator)

include_directories(2pi_event_generator)
include_directories(2pi_event_generator/get_xsect)
include_directories(2pi_event_generator/interpol)
add_executable(twopeg.exe 2pi_event_generator/main_prog.cxx)
target_link_libraries(twopeg.exe ${ROOT_LIBRARIES} 2pilib interpol getXsect)
