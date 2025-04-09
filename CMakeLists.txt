cmake_minimum_required(VERSION 3.28)
project(EPSMIM)

set(CMAKE_CXX_STANDARD 23)

#SET(CMAKE_CXX_COMPILER clang++)

find_package(likwid REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(src/affinity src/realtime_plot)

set(COMMON_SRC src/realtime_plot/GnuplotRealTime.cpp src/affinity/affinity.cpp)

add_executable(v1 src/v1/main.cpp
        src/v1/wave_solver.cpp
        ${COMMON_SRC})

add_executable(v1.1 src/v1.1/main.cpp
        src/v1.1/wave_solver.cpp
        ${COMMON_SRC})

add_executable(v1.2 src/v1.2/main.cpp
        src/v1.2/wave_solver.cpp
        ${COMMON_SRC})

add_executable(v1.3 src/v1.3/main.cpp
        src/v1.3/wave_solver.cpp
        ${COMMON_SRC})

add_executable(v1.4 src/v1.4/main.cpp
        src/v1.4/wave_solver.cpp
        ${COMMON_SRC})

add_executable(v1.5 src/v1.5/main.cpp
        src/v1.5/wave_solver.cpp
        ${COMMON_SRC})

# auto compiler simd
add_executable(v2 src/v2/main.cpp
        src/v2/wave_solver.cpp
        ${COMMON_SRC})

# __restrict__ hint
add_executable(v2.1 src/v2.1/main.cpp
        src/v2.1/wave_solver.cpp
        ${COMMON_SRC})

# unaligned intrinsics without fma
add_executable(v2.2 src/v2.2/main.cpp
        src/v2.2/wave_solver.cpp
        ${COMMON_SRC})

# unaligned intrinsics with fma
add_executable(v2.3 src/v2.3/main.cpp
        src/v2.3/wave_solver.cpp
        ${COMMON_SRC})

# assume aligned without fma
add_executable(v2.4 src/v2.4/main.cpp
        src/v2.4/wave_solver.cpp
        ${COMMON_SRC}
        src/v2.4/align_alloc.h)

target_compile_options(v1 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.1 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.2 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.3 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.4 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v1.5 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v2 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v2.1 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v2.2 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)
target_compile_options(v2.3 PRIVATE --warn-no-missing-field-initializers -Wextra -Wshadow -Wpedantic)

target_include_directories(v1.5 PRIVATE ${LIKWID_INCLUDE_DIRS}) # these are for likwid benchmarking
target_link_libraries(v1.5 PRIVATE ${LIKWID_LIBRARIES})
#target_include_directories(v1.5 PRIVATE ${LIKWID_INCLUDE_DIRS})
#target_link_libraries(v1.5 PRIVATE ${LIKWID_LIBRARIES})

target_compile_definitions(v1.5 PRIVATE LIKWID_PERFMON)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -g -DNDEBUG -march=native -fopt-info-vec")
#no vector test
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -g -DNDEBUG -march=native -fno-tree-vectorize -fopt-info-vec")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -g -DNDEBUG -march=native -fopt-info-vec -fno-strict-aliasing")

message(STATUS "Compile flags for release: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "Compile flags for debug: ${CMAKE_CXX_FLAGS_DEBUG}")
