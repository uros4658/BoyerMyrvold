cmake_minimum_required(VERSION 3.15)
project(Implementation)

set(CMAKE_CXX_STANDARD 17)

add_executable(Implementation
        src/main.cpp
        src/BoyerMyrvold.cpp
        src/Graph.cpp
        src/BoyerMyrvold.cpp
        include/FaceHandle.h
        src/FaceHandle.cpp
        include/Policies.h
        include/face_vertex_iterator.h
)