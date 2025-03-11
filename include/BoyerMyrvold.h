#ifndef BOYERMYRVOLD_H
#define BOYERMYRVOLD_H

#include "Graph.h"
#include <vector>
#include <stack>
#include <map>

class BoyerMyrvold {
public:
    BoyerMyrvold(Graph& graph);
    bool isPlanar();
    bool isK5();
    bool hasK33Subgraph();
    int edgeCount();
    void findBiconnectedComponents();

private:
    Graph& g;
    std::vector<bool> visited;
    std::vector<int> dfsNum;
    std::vector<int> dfsLow;
    std::vector<int> parent;
    std::vector<std::vector<int>> biconnectedComponents;
    std::stack<std::pair<int, int>> edgeStack;
    int time;

    bool isPlanarComponent(const Graph& component);
    bool tryPlanarEmbedding(const Graph& component);
    void dfs(const Graph& component, int u, std::vector<bool>& visited, std::vector<int>& dfsOrder, std::vector<int>& parent);
    void addEdgeToEmbedding(std::vector<std::vector<int>>& embedding, int u, int v);
    bool canAddEdgeToEmbedding(std::vector<std::vector<int>>& embedding, int u, int v);
    bool doEdgesIntersect(int u1, int v1, int u2, int v2);
    bool doSegmentsIntersect(std::pair<int, int> p1, std::pair<int, int> q1, std::pair<int, int> p2, std::pair<int, int> q2);
    bool onSegment(std::pair<int, int> p, std::pair<int, int> q, std::pair<int, int> r);
    std::vector<std::vector<int>> findFaces(const std::vector<std::vector<int>>& embedding);
    void dfsForBiconnected(int u, std::stack<std::pair<int, int>>& edgeStack);
    bool checkForK33Structure(const Graph& subgraph);

    std::map<int, std::pair<int, int>> vertexCoordinates; // Assumed to be defined elsewhere
};

#endif // BOYERMYRVOLD_H