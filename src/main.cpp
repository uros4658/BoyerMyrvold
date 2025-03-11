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

    if (!isPlanar) {
        cout << "Contains K5 subgraph: " << (bm.isK5() ? "Yes" : "No") << endl;
        cout << "Contains K3,3 subgraph: " << (bm.hasK33Subgraph() ? "Yes" : "No") << endl;
    } else {
        cout << "Finding biconnected components..." << endl;
        BiconnectedComponents bc(g);
        bc.findComponents();
        cout << "Number of biconnected components: " << bc.components.size() << endl;
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
    // Connect each vertex from first part (0,1,2) to each vertex in second part (3,4,5)
    for (int i = 0; i < 3; i++) {
        for (int j = 3; j < 6; j++) {
            g.addEdge(i, j);
        }
    }
    return g;
}

Graph createPlanarGraph() {
    Graph g(8);
    // Create a cube graph (which is planar)
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
    // Test with K5 (non-planar)
    Graph k5 = createK5();
    testPlanarity("K5 (Complete graph with 5 vertices)", k5);

    // Test with K3,3 (non-planar)
    Graph k33 = createK33();
    testPlanarity("K3,3 (Complete bipartite graph)", k33);

    // Test with a planar graph
    Graph planar = createPlanarGraph();
    testPlanarity("Cube graph (Planar)", planar);

    return 0;
}

// TODO:
// 1. Implement the full Boyer-Myrvold algorithm with:
//    - Proper embedding construction that returns the actual planar embedding
//    - Kuratowski subgraph extraction for non-planar graphs (K5 or K3,3 subdivision)
//    - DFS-based vertex addition in the correct order
//    - Left-Right planarity test approach for biconnected components
//    - Walkup and Walkdown procedures for efficient path finding
//
// 2. Improve the current implementation:
//    - Replace geometric edge intersection with topological embedding
//    - Remove dependency on vertex coordinates which aren't part of the algorithm
//    - Implement proper face tracking during embedding
//
// 3. Add visualization features:
//    - Output the planar embedding in a visual format
//    - Visualize the Kuratowski subgraph when non-planar
//
// 4. Add test cases:
//    - More complex planar and non-planar graphs
//    - Performance testing with large graphs
//    - Edge cases like disconnected graphs
//
// 5. Optimization:
//    - Ensure linear-time complexity O(n) as per the original algorithm
//    - Memory optimization for large graphs