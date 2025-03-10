#ifndef GRAPH_H
#define GRAPH_H

#include <vector>

using namespace std;

class Graph {
public:
    int V;  // Number of vertices
    vector<vector<int>> adj;  // Adjacency list

    Graph(int V);
    void addEdge(int u, int v);
};

#endif