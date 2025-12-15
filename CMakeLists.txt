cmake_minimum_required(VERSION 3.18)
project(LBMtest LANGUAGES CXX)

# ---- C++ standard ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- Kokkos ----
find_package(Kokkos REQUIRED)

# ---- Executable ----
add_executable(LBMtest
    main.cpp
    ldc3dt.cpp
)

target_include_directories(LBMtest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# ---- Compile definitions (from settings.hpp) ----
target_compile_definitions(LBMtest PRIVATE
    USE_DOUBLE
    SOA
)

# ---- Link libraries ----
target_link_libraries(LBMtest PRIVATE
    Kokkos::kokkos
)
