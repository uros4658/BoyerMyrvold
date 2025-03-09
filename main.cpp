#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>

using namespace std;

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
                    while (!edgeStack.empty() && edgeStack.top() != make_pair(u, v)) {
                        component.push_back(edgeStack.top());
                        edgeStack.pop();
                    }
                    component.push_back(edgeStack.top());
                    edgeStack.pop();
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
            }
        }
    }
};

// ----------------------------- Face Merging (Union-Find) -----------------------------

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
            parent[u] = find(parent[u]);
        }
        return parent[u];
    }

    void merge(int u, int v) {
        int rootU = find(u);
        int rootV = find(v);
        if (rootU != rootV) {
            if (rank[rootU] > rank[rootV]) {
                parent[rootV] = rootU;
            } else if (rank[rootU] < rank[rootV]) {
                parent[rootU] = rootV;
            } else {
                parent[rootV] = rootU;
                rank[rootU]++;
            }
        }
    }
};

// ----------------------------- Boyer-Myrvold Planarity Test -----------------------------

class BoyerMyrvold {
public:
    Graph& g;
    vector<int> num, low, parent;
    stack<int> pathStack;
    int time;
    FaceUnionFind faceUnion;

    BoyerMyrvold(Graph& graph) : g(graph), time(0), faceUnion(graph.V) {
        num.assign(g.V, -1);
        low.assign(g.V, -1);
        parent.assign(g.V, -1);
    }

    void dfs(int u) {
        num[u] = low[u] = time++;
        pathStack.push(u);

        for (int v : g.adj[u]) {
            if (num[v] == -1) {
                parent[v] = u;
                dfs(v);
                low[u] = min(low[u], low[v]);
            } else if (v != parent[u]) {
                low[u] = min(low[u], num[v]);

                // Merge faces correctly
                if (faceUnion.find(u) == faceUnion.find(v)) {
                    throw "Graph is not planar!";
                }
                faceUnion.merge(u, v);
            }
        }
    }

    bool isPlanar() {
        for (int i = 0; i < g.V; i++) {
            if (num[i] == -1) {
                try {
                    dfs(i);
                } catch (...) {
                    return false;
                }
            }
        }
        return true;
    }
};

// ----------------------------- Main Program -----------------------------

int main() {
    Graph g(5);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(3, 4);
    g.addEdge(4, 0);

    BiconnectedComponents bc(g);
    bc.findComponents();

    for (auto& component : bc.components) {
        Graph subgraph(g.V);
        for (auto& edge : component) {
            subgraph.addEdge(edge.first, edge.second);
        }

        BoyerMyrvold bm(subgraph);
        if (!bm.isPlanar()) {
            cout << "Graph is NOT planar!" << endl;
            return 0;
        }
    }

    cout << "Graph is planar!" << endl;
    return 0;
}




// TODO: Implement a proper embedding structure (half-edge or doubly linked face representation)

// TODO: Maintain a doubly linked list of faces for incremental embedding

// TODO: Modify DFS to correctly track edge insertions into faces

// TODO: Implement path merging logic for handling back edges properly

// TODO: Ensure correct handling of lowpoints and face merging to maintain planarity

// TODO: Optimize face union-find operations for efficiency
