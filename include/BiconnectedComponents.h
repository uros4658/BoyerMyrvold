#ifndef BICONNECTED_COMPONENTS_H
#define BICONNECTED_COMPONENTS_H

#include "Graph.h"
#include <vector>
#include <stack>


class BiconnectedComponents {
public:
    Graph& graph;
    std::vector<int> disc, low, parent;
    std::stack<std::pair<int, int>> edgeStack;
    std::vector<std::vector<std::pair<int, int>>> components;
    int time;

    BiconnectedComponents(Graph& g);
    void dfs(int u);
    void findComponents();
};

#endif