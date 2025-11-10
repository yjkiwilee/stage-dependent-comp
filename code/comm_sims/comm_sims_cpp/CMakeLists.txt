cmake_minimum_required(VERSION 3.10.0)
project(comm_sim VERSION 0.1.0 LANGUAGES C CXX)

add_executable(
    comm_sim
    main.cpp
)

target_compile_features(comm_sim PUBLIC cxx_std_17)

include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
