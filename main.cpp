#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>
#include <stdexcept>
#include <functional>
using namespace std;

// Graph class
class Graph {
public:
    int V;
    vector<vector<int>> adj;

    Graph(int V) : V(V) {
        adj.resize(V);
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
};

// ----------------------------- Biconnected Components (Tarjan's Algorithm) -----------------------------
class BiconnectedComponents {
public:
    Graph& graph;
    vector<int> disc, low, parent;
    stack<pair<int, int>> edgeStack;
    vector<vector<pair<int, int>>> components;
    int time;

    BiconnectedComponents(Graph& g) : graph(g) {
        disc.assign(graph.V, -1);
        low.assign(graph.V, -1);
        parent.assign(graph.V, -1);
        time = 0;
    }

    void dfs(int u) {
        disc[u] = low[u] = time++;
        int children = 0;
        for (int v : graph.adj[u]) {
            if (disc[v] == -1) {
                parent[v] = u;
                edgeStack.push({u, v});
                children++;
                dfs(v);
                low[u] = min(low[u], low[v]);
                if ((parent[u] == -1 && children > 1) || (parent[u] != -1 && low[v] >= disc[u])) {
                    vector<pair<int, int>> component;
                    while (!edgeStack.empty()) {
                        auto edge = edgeStack.top();
                        edgeStack.pop();
                        component.push_back(edge);
                        if (edge.first == u && edge.second == v) {
                            break;
                        }
                    }
                    components.push_back(component);
                }
            } else if (v != parent[u] && disc[v] < disc[u]) {
                low[u] = min(low[u], disc[v]);
                edgeStack.push({u, v});
            }
        }
    }

    void findComponents() {
        for (int i = 0; i < graph.V; i++) {
            if (disc[i] == -1) {
                dfs(i);

                // Handle any remaining edges in the stack (for the last component)
                if (!edgeStack.empty()) {
                    vector<pair<int, int>> component;
                    while (!edgeStack.empty()) {
                        component.push_back(edgeStack.top());
                        edgeStack.pop();
                    }
                    components.push_back(component);
                }
            }
        }
    }
};

// ----------------------------- Optimized Face Merging (Union-Find) -----------------------------
class FaceUnionFind {
public:
    vector<int> parent, rank;

    FaceUnionFind(int n) {
        parent.resize(n);
        rank.resize(n, 0);
        for (int i = 0; i < n; i++) {
            parent[i] = i;
        }
    }

    int find(int u) {
        if (parent[u] != u) {
            parent[u] = find(parent[u]); // Path compression
        }
        return parent[u];
    }

    bool merge(int u, int v) {
        int rootU = find(u);
        int rootV = find(v);
        if (rootU == rootV) {
            return false; // Already connected, which indicates a planarity violation
        }

        if (rank[rootU] > rank[rootV]) {
            parent[rootV] = rootU;
        } else if (rank[rootU] < rank[rootV]) {
            parent[rootU] = rootV;
        } else {
            parent[rootV] = rootU;
            rank[rootU]++;
        }
        return true;
    }
};

// ----------------------------- Boyer‑Myrvold Planarity Test -----------------------------
class BoyerMyrvold {
public:
    Graph& g;
    vector<bool> visited;
    vector<int> dfsNum, dfsLow, parent;
    vector<vector<int>> biconnectedComponents;
    stack<pair<int, int>> edgeStack;
    int time;

    BoyerMyrvold(Graph& graph) : g(graph), time(0) {
        visited.resize(g.V, false);
        dfsNum.resize(g.V, -1);
        dfsLow.resize(g.V, -1);
        parent.resize(g.V, -1);
    }

    // Check if graph is K₅ (complete graph on 5 vertices)
    bool isK5() {
        if (g.V < 5) return false;

        // Count vertices with degree 4 or more
        int highDegreeVerts = 0;
        for (int i = 0; i < g.V; i++) {
            if (g.adj[i].size() >= 4) highDegreeVerts++;
        }

        if (highDegreeVerts < 5) return false;

        // Identify first 5 vertices with degree 4+
        vector<int> candidates;
        for (int i = 0; i < g.V && candidates.size() < 5; i++) {
            if (g.adj[i].size() >= 4) candidates.push_back(i);
        }

        // Check if they form a complete subgraph
        for (int i = 0; i < 5; i++) {
            for (int j = i + 1; j < 5; j++) {
                int u = candidates[i];
                int v = candidates[j];
                if (find(g.adj[u].begin(), g.adj[u].end(), v) == g.adj[u].end()) {
                    return false; // Edge missing between these vertices
                }
            }
        }
        return true;
    }

    // Check if graph contains K₃,₃
    bool hasK33Subgraph() {
        // This is a simplified check - a complete implementation would
        // require more sophisticated algorithms
        if (g.V < 6) return false;

        // A preliminary check: K₃,₃ requires at least 6 vertices with degree ≥ 3
        int vertexCount = 0;
        for (int i = 0; i < g.V; i++) {
            if (g.adj[i].size() >= 3) vertexCount++;
        }
        return vertexCount >= 6 && edgeCount() >= 9; // K₃,₃ has 9 edges
    }

    int edgeCount() {
        int count = 0;
        for (int i = 0; i < g.V; i++) {
            count += g.adj[i].size();
        }
        return count / 2; // Each edge is counted twice
    }

    // Find biconnected components using Tarjan's algorithm
    void findBiconnectedComponents() {
        for (int i = 0; i < g.V; i++) {
            if (dfsNum[i] == -1) {
                dfsForBiconnected(i);
            }
        }
    }

    void dfsForBiconnected(int u) {
        visited[u] = true;
        dfsNum[u] = dfsLow[u] = time++;
        int children = 0;

        for (int v : g.adj[u]) {
            if (dfsNum[v] == -1) {
                children++;
                parent[v] = u;
                edgeStack.push({u, v});

                dfsForBiconnected(v);
                dfsLow[u] = min(dfsLow[u], dfsLow[v]);

                if ((parent[u] == -1 && children > 1) ||
                    (parent[u] != -1 && dfsLow[v] >= dfsNum[u])) {
                    // Found an articulation point, extract the biconnected component
                    vector<int> component;
                    pair<int, int> edge;
                    do {
                        edge = edgeStack.top();
                        edgeStack.pop();

                        // Add unique vertices to component
                        if (find(component.begin(), component.end(), edge.first) == component.end()) {
                            component.push_back(edge.first);
                        }
                        if (find(component.begin(), component.end(), edge.second) == component.end()) {
                            component.push_back(edge.second);
                        }
                    } while (edge.first != u || edge.second != v);

                    biconnectedComponents.push_back(component);
                }
            }
            else if (v != parent[u] && dfsNum[v] < dfsNum[u]) {
                // Back edge
                edgeStack.push({u, v});
                dfsLow[u] = min(dfsLow[u], dfsNum[v]);
            }
        }
    }

    bool isPlanar() {
        // Direct test for K₅
        if (isK5()) {
            return false;
        }

        // Direct test for K₃,₃
        if (hasK33Subgraph()) {
            return false;
        }

        // Find biconnected components to test each separately
        findBiconnectedComponents();

        // Apply Euler's formula check for each biconnected component
        for (auto& component : biconnectedComponents) {
            // Extract the subgraph for this component
            unordered_map<int, int> vertexMap;
            for (int i = 0; i < component.size(); i++) {
                vertexMap[component[i]] = i;
            }

            Graph subgraph(component.size());
            for (int u : component) {
                for (int v : g.adj[u]) {
                    if (find(component.begin(), component.end(), v) != component.end()) {
                        subgraph.addEdge(vertexMap[u], vertexMap[v]);
                    }
                }
            }

            // Check Euler's formula: E <= 3V - 6 for a maximal planar graph
            int V = subgraph.V;
            int E = 0;
            for (int i = 0; i < V; i++) {
                E += subgraph.adj[i].size();
            }
            E /= 2; // Each edge is counted twice

            if (V >= 3 && E > 3*V - 6) {
                return false; // Too many edges for a planar graph
            }
        }

        return true; // All tests passed, graph is planar
    }
};
// ----------------------------- Embedding Structure (Half-Edge Representation) -----------------------------
class HalfEdge {
public:
    int origin;
    HalfEdge* twin;
    HalfEdge* next;
    int faceID;

    HalfEdge(int orig) : origin(orig), twin(nullptr), next(nullptr), faceID(-1) {}
};

class Face {
public:
    vector<HalfEdge*> boundary;
    Face* prev;
    Face* next;

    Face() : prev(nullptr), next(nullptr) {}

    void addEdge(HalfEdge* edge) {
        boundary.push_back(edge);
        edge->faceID = 0; // Placeholder for face ID assignment.
    }
};

class Embedding {
public:
    vector<HalfEdge*> halfEdges;
    vector<Face*> faces;
    Face* faceHead;

    // For DFS lowpoint computation.
    vector<int> disc;
    vector<int> low;
    int dfsTime;

    Embedding() : faceHead(nullptr), dfsTime(0) {}

    void addEdge(int u, int v) {
        HalfEdge* edge1 = new HalfEdge(u);
        HalfEdge* edge2 = new HalfEdge(v);
        edge1->twin = edge2;
        edge2->twin = edge1;
        halfEdges.push_back(edge1);
        halfEdges.push_back(edge2);
    }

    void addFace(Face* f) {
        if (!faceHead) {
            faceHead = f;
        } else {
            f->next = faceHead;
            faceHead->prev = f;
            faceHead = f;
        }
        faces.push_back(f);
    }

    void constructFaces() {
        // Mark all half-edges as unvisited
        vector<bool> visited(halfEdges.size(), false);

        // Create faces by following next pointers in CCW order
        for (size_t i = 0; i < halfEdges.size(); i++) {
            if (!visited[i]) {
                HalfEdge* start = halfEdges[i];
                if (!start || !start->next) continue;

                Face* face = new Face();
                HalfEdge* curr = start;
                do {
                    face->addEdge(curr);
                    // Mark this half-edge as visited
                    visited[distance(halfEdges.begin(),
                           find(halfEdges.begin(), halfEdges.end(), curr))] = true;
                    curr = curr->next;
                } while (curr && curr != start);

                addFace(face);
            }
        }
    }
    void dfsEmbeddingLowpoint(int u, vector<bool>& visited, vector<HalfEdge*>& lastEdge,
                             const vector<vector<HalfEdge*>>& outgoing) {
        visited[u] = true;
        disc[u] = dfsTime;
        low[u] = dfsTime;
        dfsTime++;

        for (HalfEdge* edge : outgoing[u]) {
            int v = edge->twin->origin;

            // Link this edge to the previous one in the embedding
            if (lastEdge[u] != nullptr) {
                lastEdge[u]->next = edge;
            }
            lastEdge[u] = edge;

            if (!visited[v]) {
                // Tree edge
                dfsEmbeddingLowpoint(v, visited, lastEdge, outgoing);
                low[u] = min(low[u], low[v]);

                // Update next pointers based on lowpoint values
                if (low[v] < disc[u] && lastEdge[u] && lastEdge[u]->next == nullptr) {
                    lastEdge[u]->next = edge;
                }
            } else {
                // Back edge
                low[u] = min(low[u], disc[v]);

                // Update next pointers for back edges
                if (lastEdge[u] && lastEdge[u]->next == nullptr) {
                    lastEdge[u]->next = edge;
                }

                if (lastEdge[v] && edge->next == nullptr) {
                    edge->next = lastEdge[v];
                }
            }
        }
    }

    void buildOutgoingMapLowpoint(int numVertices) {
        vector<vector<HalfEdge*>> outgoing(numVertices);
        for (HalfEdge* edge : halfEdges) {
            outgoing[edge->origin].push_back(edge);
        }

        disc.assign(numVertices, -1);
        low.assign(numVertices, -1);
        vector<bool> visited(numVertices, false);
        vector<HalfEdge*> lastEdge(numVertices, nullptr);

        // Run DFS on each connected component
        for (int i = 0; i < numVertices; i++) {
            if (!visited[i]) {
                dfsEmbeddingLowpoint(i, visited, lastEdge, outgoing);
            }
        }
    }

    ~Embedding() {
        // Clean up memory
        for (HalfEdge* edge : halfEdges) {
            delete edge;
        }
        for (Face* face : faces) {
            delete face;
        }
    }
};

// ----------------------------- Main Program -----------------------------
int main() {
    // Create a nonplanar graph (complete graph K_5)
    Graph g(5);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(0, 3);
    g.addEdge(0, 4);
    g.addEdge(1, 2);
    g.addEdge(1, 3);
    g.addEdge(1, 4);
    g.addEdge(2, 3);
    g.addEdge(2, 4);
    g.addEdge(3, 4);

    // Process biconnected components
    BiconnectedComponents bc(g);
    bc.findComponents();
    cout << "Found " << bc.components.size() << " biconnected components" << endl;

    // Check planarity using improved Boyer-Myrvold implementation
    BoyerMyrvold bm(g);
    if (bm.isPlanar()) {
        cout << "Graph is planar!" << endl;
    } else {
        cout << "Graph is NOT planar!" << endl;
    }

    // Embedding structure demonstration with a simple triangle
    Embedding embedding;
    embedding.addEdge(0, 1);
    embedding.addEdge(1, 2);
    embedding.addEdge(2, 0);

    // Correctly set up the next pointers to form the faces
    if (embedding.halfEdges.size() >= 6) {
        // First face (CCW orientation)
        embedding.halfEdges[0]->next = embedding.halfEdges[2];
        embedding.halfEdges[2]->next = embedding.halfEdges[4];
        embedding.halfEdges[4]->next = embedding.halfEdges[0];

        // Second face (CCW orientation)
        embedding.halfEdges[1]->next = embedding.halfEdges[5];
        embedding.halfEdges[5]->next = embedding.halfEdges[3];
        embedding.halfEdges[3]->next = embedding.halfEdges[1];
    }

    embedding.constructFaces();
    embedding.buildOutgoingMapLowpoint(3); // Triangle has 3 vertices

    cout << "Number of faces detected: " << embedding.faces.size() << endl;

    // Traverse the doubly linked list of faces
    cout << "Traversing the doubly linked list of faces:" << endl;
    Face* current = embedding.faceHead;
    int faceCount = 0;
    while (current) {
        cout << "Face " << faceCount++ << " with " << current->boundary.size() << " edges." << endl;
        current = current->next;
    }

    return 0;
}