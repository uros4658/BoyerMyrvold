#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <filesystem>
#include "../include/Graph.h"
#include "../include/BoyerMyrvold.h"
#include "../include/Policies.h"
#include "../include/FaceHandle.h"

using namespace std;
namespace fs = filesystem;

// Function to decode a g6 string into a graph
Graph decodeG6(const string& g6str) {
    if (g6str.empty()) {
        return Graph(0);
    }

    // Parse number of vertices from the first character
    int n = g6str[0] - 63;
    Graph g(n);

    // If there are no edges, return the empty graph
    if (g6str.length() == 1) {
        return g;
    }

    // Decode the edge information
    int bitPos = 0;
    for (int i = 1; i < g6str.length(); i++) {
        int byte = g6str[i] - 63;

        for (int j = 0; j < 6; j++) {
            if ((byte & (1 << (5 - j))) != 0) {
                int k = bitPos + j;
                // Convert edge number to pair of vertices
                int v = 0;
                while (k >= v) {
                    k -= v;
                    v++;
                }
                int w = k;
                g.addEdge(v, w);
            }
        }
        bitPos += 6;
    }

    return g;
}

void testPlanarity(const string& graphName, Graph& g) {
    cout << "Testing planarity of " << graphName << "..." << endl;
    cout << "Graph has " << g.V << " vertices" << endl;

    int edges = 0;
    for (int i = 0; i < g.V; i++) {
        edges += g.adj[i].size();
    }
    cout << "Graph has " << edges/2 << " edges" << endl;

    BoyerMyrvold bm(g);

    try {
        bool isPlanar = bm.isPlanar();

        if (isPlanar) {
            cout << "Graph is planar!" << endl;

            // Find and print the external face for visualization
            if (!g.adj.empty()) {
                vector<int> externalFace = bm.findExternalFace(0);
                cout << "External face vertices: ";
                for (int v : externalFace) {
                    cout << v << " ";
                }
                cout << endl;
            }
        } else {
            cout << "Graph is NOT planar!" << endl;
        }
    } catch (const exception& e) {
        cout << "Error during planarity testing: " << e.what() << endl;
    }

    cout << "----------------------------------------" << endl;
}

void testG6File(const string& filepath) {
    ifstream file(filepath);
    if (!file.is_open()) {
        cerr << "Failed to open file: " << filepath << endl;
        return;
    }

    string line;
    int count = 0;
    while (getline(file, line)) {
        if (!line.empty()) {
            cout << "Graph #" << ++count << " from " << filepath << endl;
            Graph g = decodeG6(line);
            testPlanarity("Graph from G6", g);
        }
    }

    file.close();
}

void testG6Directory(const string& dirPath) {
    try {
        for (const auto& entry : fs::directory_iterator(dirPath)) {
            if (entry.path().extension() == ".g6") {
                cout << "Testing file: " << entry.path().filename() << endl;
                testG6File(entry.path().string());
            }
        }
    } catch (const fs::filesystem_error& e) {
        cerr << "Filesystem error: " << e.what() << endl;
    }
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
    // Test predefined graphs
    cout << "Testing predefined graphs:" << endl;
    Graph k5 = createK5();
    testPlanarity("K5 (Complete graph with 5 vertices)", k5);

    Graph k33 = createK33();
    testPlanarity("K3,3 (Complete bipartite graph)", k33);

    Graph planar = createPlanarGraph();
    testPlanarity("Cube graph (Planar)", planar);

    // Test graphs from g6 files in the testGraphs directory
    cout << "\nTesting graphs from g6 files:" << endl;
    testG6Directory("../testGraphs");

    return 0;
}