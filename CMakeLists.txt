cmake_minimum_required(VERSION 3.10)
project(hestonModel)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
    include
    /opt/homebrew/include
)

# Link directories for FFTW
link_directories(
    /opt/homebrew/lib
)


add_executable(HestonModel
    src/hestonFFT.cpp
    src/sim_engine.cpp
)

target_link_libraries(HestonModel fftw3)