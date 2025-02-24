cmake_minimum_required(VERSION 3.30)
project(EPSMIM)

set(CMAKE_CXX_STANDARD 23)

add_executable(EPSMIM main.cpp
        wave_solver.cpp
        wave_solver.h
        GnuplotRealTime.cpp
        GnuplotRealTime.h)

target_compile_options(EPSMIM PRIVATE -Ofast)
