#include "../include/BoyerMyrvold.h"
#include <algorithm>
#include <map>
#include <set>
#include <stack>
using namespace std;

BoyerMyrvold::BoyerMyrvold(Graph& graph) : g(graph) {
    visited.resize(g.V, false);
    dfsNum.resize(g.V, 0);
    dfsLow.resize(g.V, 0);
    parent.resize(g.V, -1);
    time = 0;
}

bool BoyerMyrvold::isK5() {
    // A K5 subgraph requires at least 5 vertices
    if (g.V < 5) return false;

    // Use a more efficient approach to find K5
    // K5 requires every vertex to have degree at least 4
    vector<int> highDegVertices;
    for (int i = 0; i < g.V; i++) {
        if (g.adj[i].size() >= 4) {
            highDegVertices.push_back(i);
            if (highDegVertices.size() >= 5) {
                // We have 5+ candidates, check if they form K5
                break;
            }
        }
    }

    // Not enough high-degree vertices
    if (highDegVertices.size() < 5) return false;

    // Try to find a K5 using the candidates
    for (size_t i1 = 0; i1 < highDegVertices.size(); i1++) {
        for (size_t i2 = i1+1; i2 < highDegVertices.size(); i2++) {
            for (size_t i3 = i2+1; i3 < highDegVertices.size(); i3++) {
                for (size_t i4 = i3+1; i4 < highDegVertices.size(); i4++) {
                    for (size_t i5 = i4+1; i5 < highDegVertices.size(); i5++) {
                        int v1 = highDegVertices[i1];
                        int v2 = highDegVertices[i2];
                        int v3 = highDegVertices[i3];
                        int v4 = highDegVertices[i4];
                        int v5 = highDegVertices[i5];

                        // Check if these 5 vertices form a complete graph
                        if (areConnected(v1, v2) && areConnected(v1, v3) && areConnected(v1, v4) && areConnected(v1, v5) &&
                            areConnected(v2, v3) && areConnected(v2, v4) && areConnected(v2, v5) &&
                            areConnected(v3, v4) && areConnected(v3, v5) &&
                            areConnected(v4, v5)) {
                            return true;
                        }
                    }
                }
            }
        }
    }

    return false;
}

// Helper method to check if two vertices are connected
bool BoyerMyrvold::areConnected(int u, int v) {
    for (int w : g.adj[u]) {
        if (w == v) return true;
    }
    return false;
}

bool BoyerMyrvold::hasK33Subgraph() {
    // K3,3 requires at least 6 vertices
    if (g.V < 6) return false;

    // A vertex in K3,3 has degree at least 3
    vector<int> candidates;
    for (int i = 0; i < g.V; i++) {
        if (g.adj[i].size() >= 3) {
            candidates.push_back(i);
        }
    }

    // Need at least 6 candidates
    if (candidates.size() < 6) return false;

    // Try all possible ways to partition the candidates into two sets of 3
    int n = candidates.size();
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                // First set: candidates[i], candidates[j], candidates[k]
                vector<int> setA = {candidates[i], candidates[j], candidates[k]};

                // Find all vertices connected to all three vertices in setA
                vector<int> potentialSetB;
                for (int v = 0; v < g.V; v++) {
                    if (find(setA.begin(), setA.end(), v) != setA.end())
                        continue; // Skip vertices in setA

                    bool connectedToAll = true;
                    for (int u : setA) {
                        if (!areConnected(u, v)) {
                            connectedToAll = false;
                            break;
                        }
                    }
                    if (connectedToAll) {
                        potentialSetB.push_back(v);
                    }
                }

                // Check if we can form setB with at least 3 vertices
                if (potentialSetB.size() >= 3) {
                    // Try all combinations of 3 vertices for setB
                    for (size_t i2 = 0; i2 < potentialSetB.size(); i2++) {
                        for (size_t j2 = i2 + 1; j2 < potentialSetB.size(); j2++) {
                            for (size_t k2 = j2 + 1; k2 < potentialSetB.size(); k2++) {
                                vector<int> setB = {potentialSetB[i2], potentialSetB[j2], potentialSetB[k2]};

                                // Check if each vertex in setB is connected to all vertices in setA
                                bool isK33 = true;
                                for (int b : setB) {
                                    for (int a : setA) {
                                        if (!areConnected(a, b)) {
                                            isK33 = false;
                                            break;
                                        }
                                    }
                                    if (!isK33) break;
                                }

                                if (isK33) return true;
                            }
                        }
                    }
                }
            }
        }
    }

    return false;
}

bool BoyerMyrvold::checkForK33Structure(const Graph& subgraph) {
    // Check if the given subgraph has the K3,3 structure
    if (subgraph.V != 6) return false;

    vector<int> partA = {0, 1, 2};
    vector<int> partB = {3, 4, 5};

    // Each vertex in part A should connect to all vertices in part B
    for (int a : partA) {
        int connections = 0;
        for (int neighbor : subgraph.adj[a]) {
            if (find(partB.begin(), partB.end(), neighbor) != partB.end()) {
                connections++;
            }
        }
        if (connections != 3) return false;
    }

    // Each vertex in part B should connect to all vertices in part A
    for (int b : partB) {
        int connections = 0;
        for (int neighbor : subgraph.adj[b]) {
            if (find(partA.begin(), partA.end(), neighbor) != partA.end()) {
                connections++;
            }
        }
        if (connections != 3) return false;
    }

    return true;
}

int BoyerMyrvold::edgeCount() {
    int count = 0;
    for (int i = 0; i < g.V; i++) {
        count += g.adj[i].size();
    }
    return count / 2; // Each edge is counted twice
}

void BoyerMyrvold::findBiconnectedComponents() {
    visited.assign(g.V, false);
    dfsNum.assign(g.V, 0);
    dfsLow.assign(g.V, 0);
    parent.assign(g.V, -1);
    biconnectedComponents.clear();

    while (!edgeStack.empty()) edgeStack.pop(); // Clear existing stack
    time = 0;

    for (int i = 0; i < g.V; i++) {
        if (!visited[i]) {
            dfsForBiconnected(i, edgeStack);
        }
    }
}

void BoyerMyrvold::dfsForBiconnected(int u, stack<pair<int, int>>& edgeStack) {
    visited[u] = true;
    dfsNum[u] = dfsLow[u] = ++time;
    int children = 0;

    for (int v : g.adj[u]) {
        if (!visited[v]) {
            children++;
            parent[v] = u;
            edgeStack.push({u, v});

            dfsForBiconnected(v, edgeStack);

            dfsLow[u] = min(dfsLow[u], dfsLow[v]);

            if ((parent[u] == -1 && children > 1) ||
                (parent[u] != -1 && dfsLow[v] >= dfsNum[u])) {
                vector<int> component;
                pair<int, int> edge;
                do {
                    edge = edgeStack.top();
                    edgeStack.pop();
                    component.push_back(edge.first);
                    component.push_back(edge.second);
                } while (edge.first != u || edge.second != v);

                biconnectedComponents.push_back(component);
            }
        }
        else if (v != parent[u] && dfsNum[v] < dfsNum[u]) {
            dfsLow[u] = min(dfsLow[u], dfsNum[v]);
            edgeStack.push({u, v});
        }
    }
}

bool BoyerMyrvold::isPlanar() {
    // Kuratowski's theorem: A graph is planar if and only if it does not contain
    // a subgraph that is a subdivision of K5 or K3,3
    if (isK5() || hasK33Subgraph()) {
        return false;
    }

    // Edge count check
    if (g.V >= 3 && edgeCount() > 3*g.V - 6) {
        return false; // Too many edges for a planar graph
    }

    // Find biconnected components and check each one
    findBiconnectedComponents();
    for (const auto& component : biconnectedComponents) {
        // Create a subgraph for the component
        set<int> vertices;
        for (int v : component) {
            vertices.insert(v);
        }

        Graph subgraph(vertices.size());
        map<int, int> oldToNew;
        int idx = 0;
        for (int v : vertices) {
            oldToNew[v] = idx++;
        }

        for (size_t i = 0; i < component.size(); i += 2) {
            subgraph.addEdge(oldToNew[component[i]], oldToNew[component[i+1]]);
        }

        if (!isPlanarComponent(subgraph)) {
            return false;
        }
    }

    return true;
}

bool BoyerMyrvold::isPlanarComponent(const Graph& component) {
    return tryPlanarEmbedding(component);
}


void BoyerMyrvold::computeDfsAndLowpoints(int v, vector<int>& dfsNum, vector<int>& dfsParent,
    vector<int>& lowPoint, vector<bool>& visited, vector<vector<int>>& embedding,
    const Graph& component) {

    visited[v] = true;
    dfsNum[v] = dfsCounter++;
    lowPoint[v] = dfsNum[v];

    for (int w : component.adj[v]) {
        if (!visited[w]) {
            dfsParent[w] = v;
            computeDfsAndLowpoints(w, dfsNum, dfsParent, lowPoint, visited, embedding, component);
            lowPoint[v] = min(lowPoint[v], lowPoint[w]);

            // Initialize the embedding with tree edges
            embedding[v].push_back(w);
            embedding[w].push_back(v);
        } else if (w != dfsParent[v]) {
            lowPoint[v] = min(lowPoint[v], dfsNum[w]);
        }
    }
}

bool BoyerMyrvold::tryPlanarEmbedding(const Graph& component) {
    int n = component.V;
    vector<vector<int>> embedding(n);
    vector<int> dfsNum(n, -1);
    vector<int> dfsParent(n, -1);
    vector<int> lowPoint(n, -1);
    vector<bool> visited(n, false);
    vector<vector<int>> nesting(n);
    vector<vector<int>> separated(n);

    // Step 1: Perform DFS and compute lowpoints
    dfsCounter = 0;
    computeDfsAndLowpoints(0, dfsNum, dfsParent, lowPoint, visited, embedding, component);

    // Step 2: Process vertices in reverse DFS order
    vector<int> reverseDfs;
    for (int i = 0; i < n; i++) {
        if (dfsNum[i] != -1) {
            reverseDfs.push_back(i);
        }
    }
    sort(reverseDfs.begin(), reverseDfs.end(), [&](int a, int b) {
        return dfsNum[a] > dfsNum[b];
    });

    // Step 3: Perform embedding and planarity testing
    for (int v : reverseDfs) {
        if (v == 0) continue; // Skip the root

        int parent = dfsParent[v];

        // Check each back edge from v
        for (int w : component.adj[v]) {
            if (dfsParent[w] != v && dfsNum[w] < dfsNum[v]) {
                // w is an ancestor of v, this is a back edge

                // Check if this back edge can be embedded without crossing
                bool canEmbed = true;
                for (int sep : separated[v]) {
                    if (dfsNum[sep] <= dfsNum[w] && dfsNum[w] < dfsNum[v]) {
                        canEmbed = false;
                        break;
                    }
                }

                if (!canEmbed) {
                    return false; // Graph is not planar
                }

                // Add the back edge to the embedding
                embedding[v].push_back(w);
                embedding[w].push_back(v);

                // Update nesting information
                nesting[parent].push_back(w);
            }
        }

        // Merge separated paths
        for (int nest : nesting[v]) {
            separated[parent].push_back(nest);
        }
    }

    return true; // Graph is planar
}

void BoyerMyrvold::dfs(const Graph& component, int u, vector<bool>& visited,
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

void BoyerMyrvold::addEdgeToEmbedding(vector<vector<int>>& embedding, int u, int v) {
    // Add edge in the clockwise order
    embedding[u].push_back(v);
    embedding[v].push_back(u);
}

bool BoyerMyrvold::canAddEdgeToEmbedding(vector<vector<int>>& embedding, int u, int v) {
    // Check if adding edge (u,v) would create a crossing
    for (int i = 0; i < embedding.size(); ++i) {
        for (int j = 0; j < embedding[i].size(); ++j) {
            int x = i;
            int y = embedding[i][j];
            if (x == u || x == v || y == u || y == v) {
                continue; // Skip edges that share a vertex with (u,v)
            }
            if (doEdgesIntersect(u, v, x, y)) {
                return false; // Found an intersection
            }
        }
    }
    return true; // No intersections found
}

bool BoyerMyrvold::doEdgesIntersect(int u1, int v1, int u2, int v2) {
    // Helper function to check if edges (u1,v1) and (u2,v2) intersect
    // This function assumes that the vertices are in a 2D plane and uses
    // a geometric approach to check for intersection
    // For simplicity, we assume vertices are points in a plane with coordinates
    // stored in a map (vertex -> (x, y))

    pair<int, int> p1 = vertexCoordinates[u1];
    pair<int, int> q1 = vertexCoordinates[v1];
    pair<int, int> p2 = vertexCoordinates[u2];
    pair<int, int> q2 = vertexCoordinates[v2];

    return doSegmentsIntersect(p1, q1, p2, q2);
}

bool BoyerMyrvold::doSegmentsIntersect(pair<int, int> p1, pair<int, int> q1, pair<int, int> p2, pair<int, int> q2) {
    // Check if line segments (p1,q1) and (p2,q2) intersect
    auto orientation = [](pair<int, int> p, pair<int, int> q, pair<int, int> r) {
        int val = (q.second - p.second) * (r.first - q.first) - (q.first - p.first) * (r.second - q.second);
        if (val == 0) return 0; // collinear
        return (val > 0) ? 1 : 2; // clock or counterclock wise
    };

    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4) return true;

    // Special cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false; // Doesn't fall in any of the above cases
}

bool BoyerMyrvold::onSegment(pair<int, int> p, pair<int, int> q, pair<int, int> r) {
    // Check if point q lies on line segment pr
    if (q.first <= max(p.first, r.first) && q.first >= min(p.first, r.first) &&
        q.second <= max(p.second, r.second) && q.second >= min(p.second, r.second)) {
        return true;
    }
    return false;
}

vector<vector<int>> BoyerMyrvold::findFaces(const vector<vector<int>>& embedding) {
    int n = embedding.size();
    vector<vector<int>> faces;

    // Keep track of visited edges to avoid traversing them twice
    vector<vector<bool>> visited(n, vector<bool>(n, false));

    for (int u = 0; u < n; u++) {
        for (int v : embedding[u]) {
            if (!visited[u][v]) {
                // Start a new face
                vector<int> face;
                int curr = u;
                int prev = -1;

                // Traverse the face
                while (true) {
                    face.push_back(curr);

                    // Mark the edge as visited in both directions
                    if (prev != -1) {
                        visited[prev][curr] = true;
                        visited[curr][prev] = true;
                    }

                    // Find the next edge in the face
                    int next = -1;
                    for (int i = 0; i < embedding[curr].size(); i++) {
                        int neighbor = embedding[curr][i];
                        if (neighbor != prev && !visited[curr][neighbor]) {
                            next = neighbor;
                            break;
                        }
                    }

                    if (next == -1 || next == u) {
                        break; // Complete the face
                    }

                    prev = curr;
                    curr = next;
                }

                // Add the face if it has at least 3 vertices
                if (face.size() >= 3) {
                    faces.push_back(face);
                }
            }
        }
    }

    return faces;
}