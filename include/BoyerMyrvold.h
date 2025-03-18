#ifndef BOYERMYRVOLD_H
#define BOYERMYRVOLD_H

#include "Graph.h"
#include <vector>
#include <stack>
#include <map>
#include <set>

class BoyerMyrvold {
public:
    BoyerMyrvold(Graph& graph);

    void initializeDataStructures();

    bool isPlanar();
    std::vector<std::vector<int>> getPlanarEmbedding();

    bool leftRightPlanarityTest(const Graph &component);

    bool isK5();
    void computeDFSOrder(const Graph &g, std::vector<int> &dfsOrder);

    bool verifyPlanarEmbedding();

    bool hasK33Subgraph();
    int edgeCount();
    void findBiconnectedComponents();
    Graph extractKuratowskiSubgraph();

private:
    Graph& g;
    std::vector<bool> visited;
    std::vector<int> dfsNum;
    std::vector<int> dfsLow;
    std::vector<int> parent;
    std::vector<std::vector<int>> biconnectedComponents;
    std::stack<std::pair<int, int>> edgeStack;
    int dfsCounter{};
    int time{};
    std::vector<std::vector<int>> embedding;

    bool findK5Subdivision(Graph& subgraph);
    bool findK33Subdivision(Graph& subgraph);
    std::vector<std::vector<int>> findPaths(int start, int end, std::vector<int>& forbidden);
    bool isPathDisjoint(const std::vector<int>& path1, const std::vector<int>& path2);
    Graph kuratowskiSubgraph;
    std::vector<bool> inEmbedding;

    void performWalkup(int v, int w, std::map<int, int>& lowpointMap, std::set<int>& pertinentRoots, std::set<int>& externallyActive);
    bool performWalkdown(int v, int root, std::set<std::pair<int, int>>& remainingEdges, const std::set<int>& externallyActive);
    std::vector<int> findExternalFace(int root);

    void computeDfsAndLowpoints(int v, std::vector<int> &dfsNum, std::vector<int> &dfsParent,
                                std::vector<int> &lowPoint,
                                std::vector<bool> &visited, std::vector<std::vector<int>> &embedding,
                                const Graph &component);

    void computeDfsWithLowpoints(const Graph &g, int u, int parent, int &counter, std::vector<int> &dfsNum,
                                 std::vector<int> &dfsParent, std::vector<int> &lowpoint1, std::vector<int> &lowpoint2,
                                 std::vector<bool> &visited, std::vector<std::vector<int>> &childList,
                                 std::vector<std::vector<int>> &backEdgeList);

    bool processBackEdge(int v, int target, const std::vector<int> &dfsNum, const std::vector<int> &dfsParent,
                         const std::vector<int> &lowpoint1, std::vector<std::vector<char>> &edgeOrientation,
                         std::vector<std::vector<int>> &pertinentRoots);

    bool mergeChildComponents(int v, int child, const std::vector<int> &dfsNum, const std::vector<int> &lowpoint1,
                              std::vector<std::vector<char>> &edgeOrientation,
                              std::vector<std::vector<int>> &pertinentRoots);

    bool addVertexToEmbedding(const Graph & graph, int i);

    bool isPlanarComponent(const Graph& component);

    bool areConnected(int u, int v);

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

    std::map<int, std::pair<int, int>> vertexCoordinates;
};

#endif // BOYERMYRVOLD_H