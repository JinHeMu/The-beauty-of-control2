cmake_minimum_required(VERSION 3.4...3.18)#设置cmake版本
project(control)#工程名字

set(CMAKE_CXX_STANDARD 17)#设置c++版本特性


include_directories("/usr/include/eigen3")

add_executable(F1_LQR_Gain #可执行程序名称
        F1_LQR_Gain.cpp
)#源文件


