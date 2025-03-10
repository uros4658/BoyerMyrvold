#ifndef BICONNECTED_COMPONENTS_H
#define BICONNECTED_COMPONENTS_H

#include "Graph.h"
#include <vector>
#include <stack>

using namespace std;

class BiconnectedComponents {
public:
    Graph& graph;
    vector<int> disc, low, parent;
    stack<pair<int, int>> edgeStack;
    vector<vector<pair<int, int>>> components;
    int time;

    BiconnectedComponents(Graph& g);
    void dfs(int u);
    void findComponents();
};

#endif