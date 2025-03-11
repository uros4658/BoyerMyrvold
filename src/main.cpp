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
// 1. Implement the `Graph` class with necessary methods like `addEdge`.
// 2. Implement the `tryPlanarEmbedding` method with a complete Boyer-Myrvold planarity test.
// 4. Implement the `findFaces` method to compute faces of the embedding.
// 5. Optimize the `isK5` and `hasK33Subgraph` methods for better performance.
// 6. Add error handling and edge cases for all methods.
// 7. Write unit tests for all methods to ensure correctness.