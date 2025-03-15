#include <iostream>
#include <vector>
#include "../include/Graph.h"
#include "../include/BiconnectedComponents.h"
#include "../include/BoyerMyrvold.h"
#include "../include/Embedding.h"
#include "../include/FaceUnionFind.h"
using namespace std;

void testPlanarity(const string& graphName, Graph& g) {
    cout << "Testing graph: " << graphName << endl;
    cout << "Number of vertices: " << g.V << endl;

    int edgeCount = 0;
    for (int i = 0; i < g.V; i++) {
        edgeCount += g.adj[i].size();
    }
    cout << "Number of edges: " << edgeCount / 2 << endl;

    BoyerMyrvold bm(g);
    bool isPlanar = bm.isPlanar();

    cout << "Is planar: " << (isPlanar ? "Yes" : "No") << endl;

    if (isPlanar) {
        auto embedding = bm.getPlanarEmbedding();
        cout << "Planar embedding:" << endl;
        for (int i = 0; i < embedding.size(); i++) {
            cout << i << ": ";
            for (int j : embedding[i]) {
                cout << j << " ";
            }
            cout << endl;
        }
    } else {
        cout << "Contains K5 subgraph: " << (bm.isK5() ? "Yes" : "No") << endl;
        cout << "Contains K3,3 subgraph: " << (bm.hasK33Subgraph() ? "Yes" : "No") << endl;
    }

    cout << "-----------------------------------" << endl;
}

Graph createK5() {
    Graph g(5);
    for (int i = 0; i < 5; i++) {
        for (int j = i + 1; j < 5; j++) {
            g.addEdge(i, j);
        }
    }
    return g;
}

Graph createK33() {
    Graph g(6);
    for (int i = 0; i < 3; i++) {
        for (int j = 3; j < 6; j++) {
            g.addEdge(i, j);
        }
    }
    return g;
}

Graph createPlanarGraph() {
    Graph g(8);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(3, 0);
    g.addEdge(4, 5);
    g.addEdge(5, 6);
    g.addEdge(6, 7);
    g.addEdge(7, 4);
    g.addEdge(0, 4);
    g.addEdge(1, 5);
    g.addEdge(2, 6);
    g.addEdge(3, 7);
    return g;
}

int main() {
    Graph k5 = createK5();
    testPlanarity("K5 (Complete graph with 5 vertices)", k5);

    Graph k33 = createK33();
    testPlanarity("K3,3 (Complete bipartite graph)", k33);

    Graph planar = createPlanarGraph();
    testPlanarity("Cube graph (Planar)", planar);

    return 0;
}

// TODO:
// 1. Complete the Boyer-Myrvold implementation:
//    - Refactor the initialization of data structures in the constructor
//    - Fix the addVertexToEmbedding, performWalkup, and performWalkdown methods
//    - Implement proper verification that the embedding is actually planar
//
// 2. Improve the current implementation:
//    - Replace geometric edge intersection with topological embedding
//    - Remove dependency on vertexCoordinates which aren't part of the algorithm
//    - Implement proper face tracking during embedding
//    - Add embedding verification based on Euler's formula (F-E+V=2)
//    - Fix biconnected components extraction for disconnected graphs
//
// 3. Add visualization features:
//    - Output the planar embedding in a visual format (GraphML, DOT, or JSON)
//    - Visualize the Kuratowski subgraph when non-planar
//    - Add circular layout drawing based on the computed embedding
//
// 4. Add test cases:
//    - More complex planar and non-planar graphs
//    - Performance testing with large graphs
//    - Edge cases like disconnected graphs
//    - Verify embedding correctness after construction
//
// 5. Optimization:
//    - Ensure linear-time complexity O(n) as per the original algorithm
//    - Memory optimization for large graphs
//    - Cache intermediate results for biconnected components
//    - Improve Kuratowski subgraph extraction algorithm