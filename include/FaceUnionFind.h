#ifndef FACE_UNION_FIND_H
#define FACE_UNION_FIND_H

#include <vector>

class FaceUnionFind {
private:
    std::vector<int> parent;
    std::vector<int> rank;

public:
    FaceUnionFind(int n);
    int find(int x);
    void unionSets(int x, int y);
    bool connected(int x, int y);
};

#endif