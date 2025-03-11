#ifndef GRAPH_H
#define GRAPH_H

#include <vector>

class Graph {
public:
    int V;  // Number of vertices
    std::vector<std::vector<int>> adj;  // Adjacency list

    Graph(int V);
    void addEdge(int u, int v);
};

#endif