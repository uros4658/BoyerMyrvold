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
    // Check if graph contains K5 (complete graph on 5 vertices)
    if (g.V < 5) return false;

    // Find a potential K5 subgraph
    for (int i = 0; i < g.V; i++) {
        if (g.adj[i].size() < 4) continue; // Not enough connections

        // Select 4 neighbors
        vector<int> potential;
        for (int neighbor : g.adj[i]) {
            potential.push_back(neighbor);
            if (potential.size() == 4) break;
        }

        if (potential.size() < 4) continue;

        // Check if these 4 neighbors form a complete graph
        bool isComplete = true;
        for (int j = 0; j < 4; j++) {
            for (int k = j+1; k < 4; k++) {
                bool connected = false;
                for (int n : g.adj[potential[j]]) {
                    if (n == potential[k]) {
                        connected = true;
                        break;
                    }
                }
                if (!connected) {
                    isComplete = false;
                    break;
                }
            }
            if (!isComplete) break;
        }

        if (isComplete) return true;
    }

    return false;
}

bool BoyerMyrvold::hasK33Subgraph() {
    // Check for K3,3 (complete bipartite graph with 3 vertices in each part)
    if (g.V < 6) return false;

    // Try to find two sets of 3 vertices each
    // This is a simplified check - would need backtracking for a complete solution
    for (int i = 0; i < g.V; i++) {
        if (g.adj[i].size() < 3) continue;

        vector<int> setA = {i};
        vector<int> setB;

        // Find potential members for set B
        for (int j : g.adj[i]) {
            setB.push_back(j);
            if (setB.size() == 3) break;
        }

        if (setB.size() < 3) continue;

        // Try to find two more vertices for set A
        for (int j = 0; j < g.V; j++) {
            if (j == i) continue;

            bool connected = true;
            for (int b : setB) {
                bool found = false;
                for (int neighbor : g.adj[j]) {
                    if (neighbor == b) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    connected = false;
                    break;
                }
            }

            if (connected) {
                setA.push_back(j);
                if (setA.size() == 3) {
                    // Check if the graph formed by setA and setB is K3,3
                    Graph subgraph(6);
                    for (int a = 0; a < 3; a++) {
                        for (int b = 0; b < 3; b++) {
                            subgraph.addEdge(a, b+3);
                        }
                    }
                    return checkForK33Structure(subgraph);
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

bool BoyerMyrvold::tryPlanarEmbedding(const Graph& component) {
    int n = component.V;
    vector<vector<int>> embedding(n);
    vector<bool> compVisited(n, false);
    vector<int> dfsOrder;
    vector<int> compParent(n, -1);

    // DFS to get ordering
    for (int i = 0; i < n; i++) {
        if (!compVisited[i]) {
            dfs(component, i, compVisited, dfsOrder, compParent);
        }
    }

    // Try adding edges in DFS order
    for (int u : dfsOrder) {
        for (int v : component.adj[u]) {
            if (canAddEdgeToEmbedding(embedding, u, v)) {
                addEdgeToEmbedding(embedding, u, v);
            } else {
                return false; // Cannot embed edge (u,v)
            }
        }
    }

    // Find faces (optional verification step)
    vector<vector<int>> faces = findFaces(embedding);

    return true;
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
    // Compute faces of the embedding
    vector<vector<int>> faces;
    // This would be a complex implementation to find all faces
    // Simplified placeholder
    return faces;
}