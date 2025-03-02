cmake_minimum_required(VERSION 3.30)
project(EPSMIM)

set(CMAKE_CXX_STANDARD 23)

add_executable(v1 src/v1/main.cpp
        src/v1/wave_solver.cpp
        src/v1/wave_solver.h
        src/realtime_plot/GnuplotRealTime.cpp
        src/realtime_plot/GnuplotRealTime.h)

add_executable(v1.1 src/v1.1/main.cpp
        src/v1.1/wave_solver.cpp
        src/v1.1/wave_solver.h
        src/realtime_plot/GnuplotRealTime.cpp
        src/realtime_plot/GnuplotRealTime.h)

add_executable(v1.2 src/v1.2/main.cpp
        src/v1.2/wave_solver.cpp
        src/v1.2/wave_solver.h
        src/realtime_plot/GnuplotRealTime.cpp
        src/realtime_plot/GnuplotRealTime.h)

add_executable(v1.3 src/v1.3/main.cpp
        src/v1.3/wave_solver.cpp
        src/v1.3/wave_solver.h
        src/realtime_plot/GnuplotRealTime.cpp
        src/realtime_plot/GnuplotRealTime.h)

target_compile_options(v1 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.1 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.2 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.3 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -g -DNDEBUG -march=native")

message(STATUS "Compile flags for release: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "Compile flags for debug: ${CMAKE_CXX_FLAGS_DEBUG}")
