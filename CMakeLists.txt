# CMAKE version of Ubuntu 16.04 is 3.5.
cmake_minimum_required(VERSION 3.5)
project(cpp_test)

if(MSVC)
  
else(MSVC)
	set(LINKLIBS
		m
    )
endif()

set(INCL 
  WAV
  cpp/
  )
set(COMPILE_OPTION)
set(SRC
  test.cpp
  cpp/HannWindow.cpp
  cpp/Ooura_FFT.cpp
  cpp/PostProcessor.cpp
  cpp/STFT.cpp
  cpp/HannWindow.h
  cpp/Ooura_FFT.h
  cpp/PostProcessor.h
  cpp/STFT.h
  WAV/WAV.h
  )

add_executable(${PROJECT_NAME} ${SRC})
target_link_libraries(${PROJECT_NAME} ${LINKLIBS})
target_include_directories(${PROJECT_NAME}	PRIVATE	${INCL}  )
target_compile_options(${PROJECT_NAME} PRIVATE ${COMPILE_OPTION} )
set_target_properties(${PROJECT_NAME} PROPERTIES  CXX_STANDARD 14 )

