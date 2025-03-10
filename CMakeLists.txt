cmake_minimum_required(VERSION 3.15)
project(Implementation)

set(CMAKE_CXX_STANDARD 17)

add_executable(Implementation
        src/main.cpp
        src/BiconnectedComponents.cpp
        src/BoyerMyrvold.cpp
        src/Embedding.cpp
        src/FaceUnionFind.cpp
        src/Graph.cpp
        src/BoyerMyrvold.cpp
)