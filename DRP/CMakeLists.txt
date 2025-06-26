# CMAKE version of Ubuntu 16.04 is 3.5.
cmake_minimum_required(VERSION 3.5)
project(DRP)

set(SRC)
set(INCL)
set(LINKLIBS)
set(COMPILE_OPTION)

include(${CMAKE_CURRENT_LIST_DIR}/../../WAV/WAV.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/DRP.cmake)


list(APPEND INCL
  ${CMAKE_CURRENT_LIST_DIR}/ )


add_executable(test ${SRC}
${CMAKE_CURRENT_LIST_DIR}/test.cpp
)


target_include_directories(test PRIVATE    ${INCL}  )
target_link_libraries(test ${LINKLIBS} )
target_compile_options(test PRIVATE ${COMPILE_OPTION} )
set_target_properties(test PROPERTIES  CXX_STANDARD 17 )

