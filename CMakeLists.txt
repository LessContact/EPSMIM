cmake_minimum_required(VERSION 3.28)
project(EPSMIM)

set(CMAKE_CXX_STANDARD 23)

include_directories(src/affinity src/realtime_plot)

add_executable(v1 src/v1/main.cpp
        src/v1/wave_solver.cpp
        src/v1/wave_solver.h
        src/realtime_plot/GnuplotRealTime.cpp
        src/realtime_plot/GnuplotRealTime.h
        src/affinity/affinity.cpp
        src/affinity/affinity.h)

add_executable(v1.1 src/v1.1/main.cpp
        src/v1.1/wave_solver.cpp
        src/v1.1/wave_solver.h
        src/realtime_plot/GnuplotRealTime.cpp
        src/realtime_plot/GnuplotRealTime.h
        src/affinity/affinity.cpp
        src/affinity/affinity.h)

add_executable(v1.2 src/v1.2/main.cpp
        src/v1.2/wave_solver.cpp
        src/v1.2/wave_solver.h
        src/realtime_plot/GnuplotRealTime.cpp
        src/realtime_plot/GnuplotRealTime.h
        src/affinity/affinity.cpp
        src/affinity/affinity.h)

add_executable(v1.3 src/v1.3/main.cpp
        src/v1.3/wave_solver.cpp
        src/v1.3/wave_solver.h
        src/realtime_plot/GnuplotRealTime.cpp
        src/realtime_plot/GnuplotRealTime.h
        src/affinity/affinity.cpp
        src/affinity/affinity.h)

add_executable(v1.4 src/v1.4/main.cpp
        src/v1.4/wave_solver.cpp
        src/v1.4/wave_solver.h
        src/realtime_plot/GnuplotRealTime.cpp
        src/realtime_plot/GnuplotRealTime.h
        src/affinity/affinity.cpp
        src/affinity/affinity.h)

add_executable(v1.5 src/v1.5/main.cpp
        src/v1.5/wave_solver.cpp
        src/v1.5/wave_solver.h
        src/realtime_plot/GnuplotRealTime.cpp
        src/realtime_plot/GnuplotRealTime.h
        src/affinity/affinity.cpp
        src/affinity/affinity.h)

target_compile_options(v1 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.1 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.2 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.3 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.4 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.5 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -g -DNDEBUG -march=native")

message(STATUS "Compile flags for release: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "Compile flags for debug: ${CMAKE_CXX_FLAGS_DEBUG}")
