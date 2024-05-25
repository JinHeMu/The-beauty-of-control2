cmake_minimum_required(VERSION 3.4...3.18)#设置cmake版本
project(control)#工程名字

set(CMAKE_CXX_STANDARD 17)#设置c++版本特性

find_package(Pangolin REQUIRED)
include_directories(${Pangolin_INCLUDE_DIRS})
include_directories("/usr/include/eigen3")

add_executable(DP_Numerical_test #可执行程序名称
        DP_Numerical_test.cpp
)#源文件

target_link_libraries(DP_Numerical_test ${Pangolin_LIBRARIES})

