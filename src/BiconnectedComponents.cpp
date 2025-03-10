#include "../include/BiconnectedComponents.h"
#include <algorithm>

using namespace std;

BiconnectedComponents::BiconnectedComponents(Graph& g) : graph(g) {
    disc.assign(graph.V, -1);
    low.assign(graph.V, -1);
    parent.assign(graph.V, -1);
    time = 0;
}

void BiconnectedComponents::dfs(int u) {
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

void BiconnectedComponents::findComponents() {
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