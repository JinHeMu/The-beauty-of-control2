cmake_minimum_required(VERSION 3.4...3.18)#设置cmake版本
project(control)#工程名字

set(CMAKE_CXX_STANDARD 17)#设置c++版本特性

set(Eigen3_DIR D:/c-package/eigen/cmake)#设置opencv路径

find_package(Eigen3 REQUIRED)



add_executable(test #可执行程序名称
    DP_Numerial_test.cpp
        )#源文件

