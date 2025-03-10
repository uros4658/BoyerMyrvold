#ifndef BOYER_MYRVOLD_H
#define BOYER_MYRVOLD_H

#include "Graph.h"
#include <vector>
#include <stack>
#include <set>

using namespace std;

class BoyerMyrvold {
public:
    Graph& g;
    vector<bool> visited;
    vector<int> dfsNum, dfsLow, parent;
    vector<vector<int>> biconnectedComponents;
    stack<pair<int, int>> edgeStack;
    int time;

    BoyerMyrvold(Graph& graph);
    bool isK5();
    bool hasK33Subgraph();
    bool checkForK33Structure(const Graph& subgraph);
    int edgeCount();
    void findBiconnectedComponents();
    bool isPlanar();
    bool isPlanarComponent(const Graph& component);
    bool tryPlanarEmbedding(const Graph& component);
    void dfsForBiconnected(int u, stack<pair<int, int>>& edgeStack);
    void dfs(const Graph& component, int u, vector<bool>& visited, vector<int>& dfsOrder, vector<int>& parent);
    void addEdgeToEmbedding(vector<vector<int>>& embedding, int u, int v);
    bool canAddEdgeToEmbedding(vector<vector<int>>& embedding, int u, int v);
    vector<vector<int>> findFaces(const vector<vector<int>>& embedding);
};

#endif