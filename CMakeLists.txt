cmake_minimum_required(VERSION 3.30)
project(EPSMIM)

set(CMAKE_CXX_STANDARD 23)

add_executable(EPSMIM main.cpp
        wave_solver.cpp
        wave_solver.h
        GnuplotRealTime.cpp
        GnuplotRealTime.h)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -finline-functions")

#target_compile_options(EPSMIM PRIVATE -Ofast)
