#include "../include/BoyerMyrvold.h"
#include <algorithm>
#include <map>
#include <set>
#include <stack>
#include <queue>
#include <functional>

using namespace std;


BoyerMyrvold::BoyerMyrvold(Graph& graph) : g(graph), kuratowskiSubgraph(0) {
    visited.resize(g.V, false);
    dfsNum.resize(g.V, 0);
    dfsLow.resize(g.V, 0);
    parent.resize(g.V, -1);
    time = 0;
    embedding.resize(g.V); // Initialize embedding
    initializeDataStructures(); // Add this line
}

void BoyerMyrvold::initializeDataStructures() {
    // Initialize vertex sets
    int n = g.V;
    inEmbedding.resize(n, false);

    // The Boyer-Myrvold algorithm doesn't need geometric coordinates
    // So we can remove the dependency on vertexCoordinates
}

Graph BoyerMyrvold::extractKuratowskiSubgraph() {
    if (!isPlanar()) {
        // First try to find a K5 subdivision
        Graph subgraph(g.V);
        if (findK5Subdivision(subgraph)) {
            return subgraph;
        }

        // If not found, try to find a K3,3 subdivision
        subgraph = Graph(g.V); // Reset subgraph
        if (findK33Subdivision(subgraph)) {
            return subgraph;
        }
    }
    return Graph(0); // Empty graph if no Kuratowski subgraph found
}

bool BoyerMyrvold::findK5Subdivision(Graph& subgraph) {
    // Try to find 5 vertices that can be connected by disjoint paths
    for (int v1 = 0; v1 < g.V; v1++) {
        if (g.adj[v1].size() < 4) continue; // Need at least degree 4

        for (int v2 : g.adj[v1]) {
            for (int v3 : g.adj[v1]) {
                if (v3 <= v2) continue;
                for (int v4 : g.adj[v1]) {
                    if (v4 <= v3) continue;
                    for (int v5 : g.adj[v1]) {
                        if (v5 <= v4) continue;

                        // Try to find disjoint paths between all pairs of these vertices
                        vector<int> vertices = {v1, v2, v3, v4, v5};
                        bool foundSubdivision = true;

                        for (int i = 0; i < 5; i++) {
                            for (int j = i + 1; j < 5; j++) {
                                // Skip path from v1 to others as we know they're adjacent
                                if (i == 0) continue;

                                vector<int> forbidden;
                                for (int k = 0; k < 5; k++) {
                                    if (k != i && k != j) forbidden.push_back(vertices[k]);
                                }

                                auto paths = findPaths(vertices[i], vertices[j], forbidden);
                                if (paths.empty()) {
                                    foundSubdivision = false;
                                    break;
                                }

                                // Add the path to the subgraph
                                for (int k = 0; k < paths[0].size() - 1; k++) {
                                    subgraph.addEdge(paths[0][k], paths[0][k + 1]);
                                }
                            }
                            if (!foundSubdivision) break;
                        }

                        if (foundSubdivision) return true;
                    }
                }
            }
        }
    }
    return false;
}

bool BoyerMyrvold::findK33Subdivision(Graph& subgraph) {
    // Try to find two sets of 3 vertices that can be connected by disjoint paths
    for (int u1 = 0; u1 < g.V; u1++) {
        for (int u2 = u1 + 1; u2 < g.V; u2++) {
            for (int u3 = u2 + 1; u3 < g.V; u3++) {
                for (int v1 = 0; v1 < g.V; v1++) {
                    if (v1 == u1 || v1 == u2 || v1 == u3) continue;

                    for (int v2 = v1 + 1; v2 < g.V; v2++) {
                        if (v2 == u1 || v2 == u2 || v2 == u3) continue;

                        for (int v3 = v2 + 1; v3 < g.V; v3++) {
                            if (v3 == u1 || v3 == u2 || v3 == u3) continue;

                            bool foundSubdivision = true;
                            vector<int> setA = {u1, u2, u3};
                            vector<int> setB = {v1, v2, v3};

                            for (int i = 0; i < 3; i++) {
                                for (int j = 0; j < 3; j++) {
                                    vector<int> forbidden;
                                    for (int k = 0; k < 3; k++) {
                                        if (k != i) forbidden.push_back(setA[k]);
                                        if (k != j) forbidden.push_back(setB[k]);
                                    }

                                    auto paths = findPaths(setA[i], setB[j], forbidden);
                                    if (paths.empty()) {
                                        foundSubdivision = false;
                                        break;
                                    }

                                    // Add the path to the subgraph
                                    for (int k = 0; k < paths[0].size() - 1; k++) {
                                        subgraph.addEdge(paths[0][k], paths[0][k + 1]);
                                    }
                                }
                                if (!foundSubdivision) break;
                            }

                            if (foundSubdivision) return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

vector<vector<int>> BoyerMyrvold::findPaths(int start, int end, vector<int>& forbidden) {
    // BFS to find paths from start to end that avoid forbidden vertices
    vector<vector<int>> result;
    vector<bool> visited(g.V, false);

    for (int v : forbidden) {
        visited[v] = true;  // Mark forbidden vertices as visited
    }

    queue<vector<int>> q;
    q.push({start});
    visited[start] = true;

    while (!q.empty() && result.size() < 3) { // Find up to 3 paths
        vector<int> path = q.front();
        q.pop();

        int current = path.back();

        if (current == end) {
            result.push_back(path);
            continue;
        }

        for (int next : g.adj[current]) {
            if (!visited[next]) {
                visited[next] = true;
                vector<int> newPath = path;
                newPath.push_back(next);
                q.push(newPath);
            }
        }
    }

    return result;
}

bool BoyerMyrvold::isPathDisjoint(const vector<int>& path1, const vector<int>& path2) {
    // Check if two paths are internally vertex-disjoint
    for (int i = 1; i < path1.size() - 1; i++) {
        for (int j = 1; j < path2.size() - 1; j++) {
            if (path1[i] == path2[j]) return false;
        }
    }
    return true;
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

// Add this method to BoyerMyrvold class
void BoyerMyrvold::computeDFSOrder(const Graph& g, vector<int>& dfsOrder) {
    int n = g.V;
    vector<bool> visited(n, false);
    dfsOrder.clear();

    // Helper function for DFS
    function<void(int)> dfs = [&](int u) {
        visited[u] = true;

        for (int v : g.adj[u]) {
            if (!visited[v]) {
                dfs(v);
            }
        }

        // Add vertex in postorder (after visiting all descendants)
        dfsOrder.push_back(u);
    };

    // Run DFS from each unvisited vertex
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            dfs(i);
        }
    }

    // dfsOrder now contains vertices in reverse DFS order
    // No need to reverse since we're adding in postorder
}

bool BoyerMyrvold::isPlanarComponent(const Graph& component) {
    // Use the Left-Right planarity test for biconnected components
    return leftRightPlanarityTest(component);
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
    if (isK5() || hasK33Subgraph()) {
        kuratowskiSubgraph = extractKuratowskiSubgraph();
        return false;
    }

    if (g.V >= 3 && edgeCount() > 3 * g.V - 6) {
        return false;
    }

    findBiconnectedComponents();
    for (const auto& component : biconnectedComponents) {
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
            subgraph.addEdge(oldToNew[component[i]], oldToNew[component[i + 1]]);
        }

        if (!isPlanarComponent(subgraph)) {
            return false;
        }
    }

    return true;
}

bool BoyerMyrvold::tryPlanarEmbedding(const Graph& component) {
    int n = component.V;
    vector<vector<int>> localEmbedding(n);
    vector<int> dfsNum(n, -1);
    vector<int> dfsParent(n, -1);
    vector<int> lowPoint(n, -1);
    vector<bool> visited(n, false);
    vector<vector<int>> nesting(n);
    vector<vector<int>> separated(n);

    dfsCounter = 0;
    computeDfsAndLowpoints(0, dfsNum, dfsParent, lowPoint, visited, localEmbedding, component);

    vector<int> reverseDfs;
    for (int i = 0; i < n; i++) {
        if (dfsNum[i] != -1) {
            reverseDfs.push_back(i);
        }
    }
    sort(reverseDfs.begin(), reverseDfs.end(), [&](int a, int b) {
        return dfsNum[a] > dfsNum[b];
    });

    for (int v : reverseDfs) {
        if (v == 0) continue;

        int parent = dfsParent[v];

        for (int w : component.adj[v]) {
            if (dfsParent[w] != v && dfsNum[w] < dfsNum[v]) {
                bool canEmbed = true;
                for (int sep : separated[v]) {
                    if (dfsNum[sep] <= dfsNum[w] && dfsNum[w] < dfsNum[v]) {
                        canEmbed = false;
                        break;
                    }
                }

                if (!canEmbed) {
                    return false;
                }

                localEmbedding[v].push_back(w);
                localEmbedding[w].push_back(v);

                nesting[parent].push_back(w);
            }
        }

        for (int nest : nesting[v]) {
            separated[parent].push_back(nest);
        }
    }

    embedding = localEmbedding; // Store the embedding
    return true;
}

vector<vector<int>> BoyerMyrvold::getPlanarEmbedding() {
    return embedding;
}

bool BoyerMyrvold::leftRightPlanarityTest(const Graph& component) {
    int n = component.V;

    // Skip trivially planar graphs
    if (n <= 4) return true;

    // Euler's formula check: |E| ≤ 3|V| - 6 for planar graphs
    int edges = 0;
    for (int i = 0; i < n; i++) {
        edges += component.adj[i].size();
    }
    edges /= 2; // Each edge is counted twice

    if (edges > 3 * n - 6) return false;

    // Initialize data structures
    vector<int> dfsNum(n, -1);
    vector<int> dfsParent(n, -1);
    vector<int> lowpoint1(n, -1);
    vector<int> lowpoint2(n, -1);
    vector<bool> visited(n, false);

    // Edge classifications
    vector<vector<int>> childList(n);
    vector<vector<int>> backEdgeList(n);

    // Edge orientation (left or right)
    vector<vector<char>> edgeOrientation(n, vector<char>(n, '?'));

    // Pertinent roots for each vertex
    vector<vector<int>> pertinentRoots(n);

    // Generate DFS ordering with lowpoints
    int dfsCounter = 0;
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            computeDfsWithLowpoints(component, i, -1, dfsCounter, dfsNum, dfsParent, lowpoint1, lowpoint2,
                visited, childList, backEdgeList);
        }
    }

    // Process vertices in reverse DFS order
    vector<int> dfsVertices(n);
    for (int i = 0; i < n; i++) {
        if (dfsNum[i] != -1) {
            dfsVertices[dfsNum[i]] = i;
        }
    }

    for (int i = n-1; i >= 0; i--) {
        if (i >= dfsVertices.size()) continue;
        int v = dfsVertices[i];

        // Process back edges from v
        for (int target : backEdgeList[v]) {
            if (!processBackEdge(v, target, dfsNum, dfsParent, lowpoint1, edgeOrientation, pertinentRoots)) {
                return false;
            }
        }

        // Merge child components
        for (int child : childList[v]) {
            if (!mergeChildComponents(v, child, dfsNum, lowpoint1, edgeOrientation, pertinentRoots)) {
                return false;
            }
        }
    }

    // Build the planar embedding based on edge orientations
    embedding.resize(n);
    for (int i = 0; i < n; i++) {
        embedding[i].clear();
        for (int j : component.adj[i]) {
            embedding[i].push_back(j);
        }

        // Sort adjacency lists based on edge orientation
        // This creates a clockwise or counterclockwise arrangement
        if (!embedding[i].empty()) {
            sort(embedding[i].begin(), embedding[i].end(), [&](int a, int b) {
                if (edgeOrientation[i][a] != edgeOrientation[i][b]) {
                    return edgeOrientation[i][a] < edgeOrientation[i][b];
                }
                return a < b;
            });
        }
    }

    return true;
}

void BoyerMyrvold::computeDfsWithLowpoints(
    const Graph& g, int u, int parent, int& counter,
    vector<int>& dfsNum, vector<int>& dfsParent,
    vector<int>& lowpoint1, vector<int>& lowpoint2,
    vector<bool>& visited,
    vector<vector<int>>& childList,
    vector<vector<int>>& backEdgeList) {

    visited[u] = true;
    dfsNum[u] = counter++;
    dfsParent[u] = parent;
    lowpoint1[u] = dfsNum[u];
    lowpoint2[u] = dfsNum[u];

    for (int v : g.adj[u]) {
        if (!visited[v]) {
            // Tree edge
            childList[u].push_back(v);

            computeDfsWithLowpoints(g, v, u, counter, dfsNum, dfsParent, lowpoint1, lowpoint2,
                                 visited, childList, backEdgeList);

            // Update lowpoints
            if (lowpoint1[v] < lowpoint1[u]) {
                lowpoint2[u] = min(lowpoint1[u], lowpoint2[v]);
                lowpoint1[u] = lowpoint1[v];
            } else if (lowpoint1[v] == lowpoint1[u]) {
                lowpoint2[u] = min(lowpoint2[u], lowpoint2[v]);
            } else {
                lowpoint2[u] = min(lowpoint2[u], lowpoint1[v]);
            }
        }
        else if (v != parent && dfsNum[v] < dfsNum[u]) {
            // Back edge
            backEdgeList[u].push_back(v);

            if (dfsNum[v] < lowpoint1[u]) {
                lowpoint2[u] = lowpoint1[u];
                lowpoint1[u] = dfsNum[v];
            } else if (dfsNum[v] != lowpoint1[u]) {
                lowpoint2[u] = min(lowpoint2[u], dfsNum[v]);
            }
        }
    }
}

bool BoyerMyrvold::processBackEdge(int v, int target, const vector<int>& dfsNum,
                                const vector<int>& dfsParent, const vector<int>& lowpoint1,
                                vector<vector<char>>& edgeOrientation,
                                vector<vector<int>>& pertinentRoots) {
    // Find lowest common ancestor (LCA)
    int lca = target;
    while (dfsNum[lca] > dfsNum[v]) {
        lca = dfsParent[lca];
    }

    // Add constraints for the path from v to LCA
    int current = v;
    while (current != lca) {
        int parent = dfsParent[current];

        // Check if we can place this edge
        if (edgeOrientation[parent][current] == 'R') {
            return false; // Conflict with existing constraint
        }

        // Mark this edge with the left orientation
        edgeOrientation[parent][current] = 'L';
        edgeOrientation[current][parent] = 'L'; // Symmetric

        // Move up the tree
        current = parent;
    }

    // Add constraints for the path from target to LCA
    current = target;
    while (current != lca) {
        int parent = dfsParent[current];

        // Check if we can place this edge
        if (edgeOrientation[parent][current] == 'L') {
            return false; // Conflict with existing constraint
        }

        // Mark this edge with the right orientation
        edgeOrientation[parent][current] = 'R';
        edgeOrientation[current][parent] = 'R'; // Symmetric

        // Move up the tree
        current = parent;
    }

    // Add target to pertinent roots for v
    pertinentRoots[v].push_back(target);

    return true;
}

bool BoyerMyrvold::mergeChildComponents(int v, int child, const vector<int>& dfsNum,
                                     const vector<int>& lowpoint1,
                                     vector<vector<char>>& edgeOrientation,
                                     vector<vector<int>>& pertinentRoots) {
    // Check for conflicts in child component orientation
    if (lowpoint1[child] < dfsNum[v]) {
        // This child has edges that go above v in the DFS tree

        // If the edge orientation is already set, verify it's consistent
        if (edgeOrientation[v][child] != '?') {
            // For consistency with lowpoint constraints, the edge must be 'L'
            if (edgeOrientation[v][child] != 'L') {
                return false; // Conflict that can't be resolved
            }
        } else {
            // Set the orientation based on lowpoint requirements
            edgeOrientation[v][child] = 'L';
            edgeOrientation[child][v] = 'L'; // Symmetric
        }
    } else {
        // The child component doesn't have edges going above v
        if (edgeOrientation[v][child] == '?') {
            // We can choose any orientation
            edgeOrientation[v][child] = 'R';
            edgeOrientation[child][v] = 'R';
        }
    }

    // Merge pertinent roots from child
    for (int root : pertinentRoots[child]) {
        pertinentRoots[v].push_back(root);
    }

    return true;
}

bool BoyerMyrvold::addVertexToEmbedding(const Graph& g, int v) {
    // Get neighbors of v that are already in the embedding (back-edges)
   vector<int> embeddedNeighbors;
    for (int u : g.adj[v]) {
        if (inEmbedding[u]) {
            embeddedNeighbors.push_back(u);
        }
    }

    // If no embedded neighbors, just add v to the embedding
    if (embeddedNeighbors.empty()) {
        inEmbedding[v] = true;
        return true;
    }

    // Initialize data structures for walkup and walkdown
    map<int, int> lowpointMap;  // Maps vertex to its lowpoint
    set<int> pertinentRoots;    // Roots of pertinent biconnected components
    set<int> externallyActive;  // Externally active vertices

    // For each embedded neighbor, perform walkup
    for (int neighbor : embeddedNeighbors) {
        performWalkup(v, neighbor, lowpointMap, pertinentRoots, externallyActive);
    }

    // Track edges that need to be embedded
    set<pair<int, int>> remainingEdges;
    for (int neighbor : embeddedNeighbors) {
        remainingEdges.insert({v, neighbor});
    }

    // For each pertinent root, perform walkdown
    for (int root : pertinentRoots) {
        if (!performWalkdown(v, root, remainingEdges, externallyActive)) {
            return false; // Not planar
        }
    }

    // Mark v as now in the embedding
    inEmbedding[v] = true;

    // All edges were successfully embedded if remainingEdges is empty
    return remainingEdges.empty();
}

void BoyerMyrvold::performWalkup(int v, int w, map<int, int>& lowpointMap, set<int>& pertinentRoots, set<int>& externallyActive) {
    int currentVertex = w;

    // Walk up the DFS tree from w until we find a pertinent root
    while (true) {
        // Get parent of current vertex
        int parentVertex = parent[currentVertex];
        if (parentVertex == -1) break; // Reached root of DFS tree

        // Update lowpoint of current vertex
        if (lowpointMap.find(currentVertex) == lowpointMap.end() ||
            lowpointMap[currentVertex] > dfsNum[w]) {
            lowpointMap[currentVertex] = dfsNum[w];
        }

        // Check if current vertex is externally active for v
        if (lowpointMap[currentVertex] < dfsNum[v]) {
            // Mark as externally active
            externallyActive.insert(currentVertex);
            // Continue walking up
            currentVertex = parentVertex;
        } else {
            // Current vertex is a pertinent root
            pertinentRoots.insert(currentVertex);
            break;
        }
    }
}

bool BoyerMyrvold::performWalkdown(int v, int root, set<pair<int, int>>& remainingEdges, const set<int>& externallyActive) {
    // Find the external face of the current embedding at the root
    vector<int> externalFace = findExternalFace(root);

    // Try to embed all remaining edges within the external face
    for (auto it = remainingEdges.begin(); it != remainingEdges.end();) {
        int neighbor = it->second;

        // Check if neighbor is on the external face
        if (find(externalFace.begin(), externalFace.end(), neighbor) != externalFace.end()) {
            // Can embed this edge
            embedding[v].push_back(neighbor);
            embedding[neighbor].push_back(v);
            it = remainingEdges.erase(it);
        } else {
            ++it;
        }
    }

    // If there are still edges to embed, the graph is not planar
    return remainingEdges.empty();
}

vector<int> BoyerMyrvold::findExternalFace(int root) {
    vector<int> face;
    set<int> visited;

    // Start with root
    face.push_back(root);
    visited.insert(root);

    // Simple face traversal (clockwise or counterclockwise)
    int current = root;
    while (true) {
        bool found = false;
        for (int neighbor : embedding[current]) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                face.push_back(neighbor);
                current = neighbor;
                found = true;
                break;
            }
        }

        if (!found || current == root) break;
    }

    return face;
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