# CMAKE version of Ubuntu 16.04 is 3.5.
cmake_minimum_required(VERSION 3.5)
project(cpp_test)

set(LINKLIBS
  m
  )
set(COMPILE_OPTION)
set(SRC 
  test.cpp 
  cpp/HannWindow.h
  cpp/Ooura_FFT.h
  cpp/PostProcessor.h
  cpp/STFT.h
  WAV/WAV.cpp
  WAV/WAV.h
  )

add_executable(${PROJECT_NAME} ${SRC})
target_link_libraries(${PROJECT_NAME} ${LINKLIBS})
target_compile_options(${PROJECT_NAME} PRIVATE ${COMPILE_OPTION} )
set_target_properties(${PROJECT_NAME} PROPERTIES  CXX_STANDARD 14 )

