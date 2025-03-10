#include <iostream>
#include <vector>
#include <stack>
#include <queue>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <set>
#include <unordered_map>
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

        // Find 5 vertices with degree at least 4
        vector<int> candidates;
        for (int i = 0; i < g.V && candidates.size() < 5; i++) {
            if (g.adj[i].size() >= 4) candidates.push_back(i);
        }

        if (candidates.size() < 5) return false;

        // Check if they form a complete subgraph
        for (int i = 0; i < 5; i++) {
            for (int j = i + 1; j < 5; j++) {
                int u = candidates[i];
                int v = candidates[j];
                // Check if v is in the adjacency list of u
                auto it = find(g.adj[u].begin(), g.adj[u].end(), v);
                if (it == g.adj[u].end()) {
                    return false; // Edge missing
                }
            }
        }
        return true;
    }
    // Improved K3,3 detection - more efficient and appropriate
    bool hasK33Subgraph() {
        // Quick check - K3,3 needs at least 6 vertices and 9 edges
        if (g.V < 6) return false;

        int totalEdges = edgeCount();
        if (totalEdges < 9) return false;

        // Check biconnected components - K3,3 must be in a single biconnected component
        findBiconnectedComponents();

        for (const auto& component : biconnectedComponents) {
            // Need at least 6 vertices in the component
            if (component.size() < 6) continue;

            // Create subgraph for this component
            unordered_map<int, int> vertexMap;
            for (int i = 0; i < component.size(); i++) {
                vertexMap[component[i]] = i;
            }

            Graph subgraph(component.size());
            int edgeCount = 0;

            // Build the subgraph
            for (int u : component) {
                for (int v : g.adj[u]) {
                    if (vertexMap.find(v) != vertexMap.end() && u < v) {
                        subgraph.addEdge(vertexMap[u], vertexMap[v]);
                        edgeCount++;
                    }
                }
            }

            // Quick edge density check (K3,3 has 9 edges for 6 vertices)
            if (edgeCount < 9) continue;

            // Check for bipartite structure with each part having 3 vertices
            if (checkForK33Structure(subgraph)) {
                return true;
            }
        }

        return false;
    }

    bool checkForK33Structure(const Graph& subgraph) {
        // Try to find a bipartite partition where each part has 3 vertices
        vector<int> partition(subgraph.V, -1); // -1 unassigned, 0 set A, 1 set B

        // Try different starting configurations
        for (int start = 0; start < subgraph.V; start++) {
            fill(partition.begin(), partition.end(), -1);

            // Start with vertex in set A
            partition[start] = 0;

            // Use BFS to build bipartite partition
            queue<int> q;
            q.push(start);

            while (!q.empty()) {
                int u = q.front();
                q.pop();

                for (int v : subgraph.adj[u]) {
                    if (partition[v] == -1) {
                        // Assign opposite partition
                        partition[v] = 1 - partition[u];
                        q.push(v);
                    } else if (partition[v] == partition[u]) {
                        // Not bipartite
                        break;
                    }
                }
            }

            // Check if we have a valid partition with 3 vertices in each set
            int countA = 0, countB = 0;
            for (int p : partition) {
                if (p == 0) countA++;
                else if (p == 1) countB++;
            }

            if (countA == 3 && countB == 3) {
                // Check if all possible edges between sets exist
                bool isK33 = true;
                for (int i = 0; i < subgraph.V; i++) {
                    if (partition[i] == 0) {
                        for (int j = 0; j < subgraph.V; j++) {
                            if (partition[j] == 1) {
                                if (find(subgraph.adj[i].begin(), subgraph.adj[i].end(), j)
                                    == subgraph.adj[i].end()) {
                                    isK33 = false;
                                    break;
                                }
                            }
                        }
                        if (!isK33) break;
                    }
                }

                if (isK33) return true;
            }
        }

        return false;
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
        biconnectedComponents.clear();
        time = 0;

        // Reset all arrays
        fill(visited.begin(), visited.end(), false);
        fill(dfsNum.begin(), dfsNum.end(), -1);
        fill(dfsLow.begin(), dfsLow.end(), -1);
        fill(parent.begin(), parent.end(), -1);

        // Empty the stack just to be safe
        while (!edgeStack.empty()) {
            edgeStack.pop();
        }

        for (int i = 0; i < g.V; i++) {
            if (dfsNum[i] == -1) {
                dfsForBiconnected(i);

                // Process any remaining edges in the stack for this connected component
                if (!edgeStack.empty()) {
                    vector<int> component;
                    while (!edgeStack.empty()) {
                        pair<int, int> edge = edgeStack.top();
                        edgeStack.pop();

                        // Add unique vertices to component
                        if (find(component.begin(), component.end(), edge.first) == component.end()) {
                            component.push_back(edge.first);
                        }
                        if (find(component.begin(), component.end(), edge.second) == component.end()) {
                            component.push_back(edge.second);
                        }
                    }

                    // Only add component if it has at least 3 vertices (minimum for a biconnected component)
                    if (component.size() >= 3) {
                        biconnectedComponents.push_back(component);
                    }
                }
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
                    } while (!edgeStack.empty() && (edge.first != u || edge.second != v));

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

        // Use Euler's formula for a quick check on the whole graph
        if (g.V >= 3) {
            int E = edgeCount();
            if (E > 3*g.V - 6) {
                return false; // Too many edges for a planar graph
            }
        }

        // Find biconnected components to test each separately
        findBiconnectedComponents();

        // Apply planarity test for each biconnected component
        for (auto& component : biconnectedComponents) {
            // Extract the subgraph for this component
            unordered_map<int, int> vertexMap;
            for (int i = 0; i < component.size(); i++) {
                vertexMap[component[i]] = i;
            }

            Graph subgraph(component.size());
            for (int u : component) {
                for (int v : g.adj[u]) {
                    if (vertexMap.find(v) != vertexMap.end() && u < v) {
                        subgraph.addEdge(vertexMap[u], vertexMap[v]);
                    }
                }
            }

            // Check if this component is planar
            if (!isPlanarComponent(subgraph)) {
                return false;
            }
        }

        return true; // All components are planar
    }

    bool isPlanarComponent(const Graph& component) {
        // For small components, planarity is guaranteed
        if (component.V <= 4) return true;

        // Perform simplified Boyer-Myrvold embedding test
        return tryPlanarEmbedding(component);
    }

    bool tryPlanarEmbedding(const Graph& component) {
        // We'll implement a simplified version of the Boyer-Myrvold algorithm
        // This performs a DFS-based embedding approach

        vector<bool> visited(component.V, false);
        vector<int> dfsOrder;
        vector<int> parentInDFS(component.V, -1);

        // Step 1: DFS traversal to establish ordering
        dfs(component, 0, visited, dfsOrder, parentInDFS);

        // Step 2: Initialize empty embedding
        vector<vector<int>> embedding(component.V);

        // Step 3: Add edges in DFS tree first
        for (int i = 0; i < component.V; i++) {
            for (int neighbor : component.adj[i]) {
                if (parentInDFS[i] == neighbor || parentInDFS[neighbor] == i) {
                    // This is a tree edge, add it to embedding
                    addEdgeToEmbedding(embedding, i, neighbor);
                }
            }
        }

        // Step 4: Add back edges while preserving planarity
        for (int i = 0; i < component.V; i++) {
            for (int neighbor : component.adj[i]) {
                if (parentInDFS[i] != neighbor && parentInDFS[neighbor] != i) {
                    // This is a back edge
                    if (!canAddEdgeToEmbedding(embedding, i, neighbor)) {
                        return false; // Cannot embed this edge without crossing
                    }
                    addEdgeToEmbedding(embedding, i, neighbor);
                }
            }
        }

        return true;
    }

    void dfs(const Graph& component, int u, vector<bool>& visited,
             vector<int>& dfsOrder, vector<int>& parent) {
        visited[u] = true;
        dfsOrder.push_back(u);

        for (int v : component.adj[u]) {
            if (!visited[v]) {
                parent[v] = u;
                dfs(component, v, visited, dfsOrder, parent);
            }
        }
    }

    void addEdgeToEmbedding(vector<vector<int>>& embedding, int u, int v) {
        // Add edge in clockwise order around each vertex
        embedding[u].push_back(v);
        embedding[v].push_back(u);
    }

    bool canAddEdgeToEmbedding(vector<vector<int>>& embedding, int u, int v) {
    // Step 1: Find all faces in the current embedding
    vector<vector<int>> faces = findFaces(embedding);

    // Step 2: Check if vertices u and v are on the same face
    for (const auto& face : faces) {
        bool containsU = false;
        bool containsV = false;

        for (int vertex : face) {
            if (vertex == u) containsU = true;
            if (vertex == v) containsV = true;
        }

        // If both vertices are on the same face, we can add the edge without crossing
        if (containsU && containsV) {
            return true;
        }
    }

    return false;
}

// Helper method to find all faces in the current embedding
vector<vector<int>> findFaces(const vector<vector<int>>& embedding) {
    vector<vector<int>> faces;

    // Keep track of visited edges to avoid duplicates
    set<pair<int, int>> visitedEdges;

    // Find faces by walking along the boundary
    for (int start = 0; start < embedding.size(); start++) {
        for (int adjVertex : embedding[start]) {
            // Skip visited edges
            if (visitedEdges.count({min(start, adjVertex), max(start, adjVertex)})) {
                continue;
            }

            // Start a new face
            vector<int> face;
            int current = start;
            int prev = -1;

            do {
                face.push_back(current);

                // Mark this edge as visited
                visitedEdges.insert({min(prev, current), max(prev, current)});

                // Find the next edge in counterclockwise order
                int next = -1;
                if (prev != -1) {
                    // Find the edge immediately counterclockwise from prev->current
                    int prevIndex = distance(embedding[current].begin(),
                                          find(embedding[current].begin(), embedding[current].end(), prev));

                    // Get the next vertex in counterclockwise order
                    next = embedding[current][(prevIndex + 1) % embedding[current].size()];
                } else {
                    // For the starting edge
                    next = adjVertex;
                }

                prev = current;
                current = next;

            } while (current != start);

            // Add the face if it's valid
            if (face.size() > 2) {
                faces.push_back(face);
            }
        }
    }

    return faces;
}
};

int main() {
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

    // Check planarity using Boyer-Myrvold implementation
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

