#include "../include/BoyerMyrvold.h"
#include "../include/FaceHandle.h"
#include "../include/Graph.h"
#include "../include/face_vertex_iterator.h"
#include "../include/planar_face_traversal_visitor.h"
#include <algorithm>
#include <set>
#include <stack>
#include <queue>
#include <functional>
#include "../include/Policies.h"

using namespace std;

BoyerMyrvold::BoyerMyrvold(Graph& graph) : g(graph) {
    initializeDataStructures();
}

StoreOldHandlesPolicy storeOldHandlesPolicy;
StoreEmbeddingPolicy storeEmbeddingPolicy;

void BoyerMyrvold::initializeDataStructures() {
    visited.resize(g.V, false);
    dfsNum.resize(g.V, 0);
    backedgeFlag.resize(g.V, -1);
    canonicalDfsChild.resize(g.V, -1);
    dfsParent.resize(g.V, -1);
    dfsChildHandles.resize(g.V);
    lowPoint.resize(g.V, -1);
    leastAncestor.resize(g.V, -1);
    pertinentRoots.resize(g.V);
    faceHandles.resize(g.V);
    flipped.resize(g.V, false);
    vertices_by_dfs_num.resize(g.V);
    backedges.resize(g.V);
    dfsParentEdge.resize(g.V);
}

bool BoyerMyrvold::isPlanar() {
    initializeDataStructures();

    // Perform DFS to initialize DFS numbering and other data structures
    int dfsCount = 0;
    vector<bool> dfsVisited(g.V, false);

    // DFS function to set up the required data
    function<void(int, int)> dfs = [&](int v, int parent) {
        dfsVisited[v] = true;
        dfsNum[v] = dfsCount++;
        dfsParent[v] = parent;
        vertices_by_dfs_num[dfsNum[v]] = v;
        lowPoint[v] = dfsNum[v];
        leastAncestor[v] = g.V; // Initialize to a high value

        // Initialize face handle for the vertex
        faceHandles[v].setFirstVertex(v);
        faceHandles[v].setSecondVertex(v);
        faceHandles[v].setAnchor(v);

        for (int w : g.adj[v]) {
            if (!dfsVisited[w]) {
                // Tree edge
                dfsParentEdge[w] = make_pair(v, w);
                dfs(w, v);
                lowPoint[v] = min(lowPoint[v], lowPoint[w]);

                // Set up DFS child handles - CRITICAL PART
                dfsChildHandles[w].setFirstVertex(v);
                dfsChildHandles[w].setSecondVertex(w);
                dfsChildHandles[w].setAnchor(w);

                // Update face handles to create initial embedding
                faceHandles[v].setFirstVertex(w);
                faceHandles[w].setFirstVertex(v);

                // Initialize canonical DFS child
                canonicalDfsChild[w] = w;
            }
            else if (w != parent) {
                // Back edge to any non-parent vertex
                if (dfsNum[w] < dfsNum[v]) {
                    // Back edge to ancestor
                    leastAncestor[v] = min(leastAncestor[v], dfsNum[w]);
                    lowPoint[v] = min(lowPoint[v], dfsNum[w]);

                    // Record backedge in both directions
                    backedges[v].push_back(make_pair(v, w));
                    backedges[w].push_back(make_pair(w, v));
                }
                // Note: Edges to later vertices are handled during walkup
            }
        }
    };

    // Start DFS from vertex 0
    for (int i = 0; i < g.V; i++) {
        if (!dfsVisited[i]) {
            dfs(i, -1);
        }
    }

    // Process vertices in reverse DFS order (bottom-up)
    for (auto vi = vertices_by_dfs_num.rbegin(); vi != vertices_by_dfs_num.rend(); ++vi) {
        int v = *vi;

        // Clear pertinent roots for current vertex before walkup
        pertinentRoots[v].clear();

        // Reset visited flags for the walkup phase
        fill(visited.begin(), visited.end(), 0);

        performWalkup(v);

        if (!performWalkdown(v)) {
            return false;
        }
    }

    return true;
}
bool BoyerMyrvold::pertinent(int vertex, int v) {
    return !backedges[vertex].empty() || backedgeFlag[vertex] == dfsNum[v];
}

bool BoyerMyrvold::externally_active(int vertex, int v) {
    for (int w : g.adj[vertex]) {
        if (dfsNum[w] < dfsNum[vertex] && w != dfsParent[vertex]) {
            return true;
        }
    }
    return false;
}

bool BoyerMyrvold::internally_active(int vertex, int v) {
    for (int w : g.adj[vertex]) {
        if (dfsNum[w] > dfsNum[vertex] && lowPoint[w] < dfsNum[vertex]) {
            return true;
        }
    }
    return false;
}

void BoyerMyrvold::remove_vertex_from_separated_dfs_child_list(int vertex) {
    for (auto it = pertinentRoots[vertex].begin(); it != pertinentRoots[vertex].end(); ++it) {
        // Use face_handle_t methods to access the first vertex
        if (it->firstVertex() == vertex) {
            pertinentRoots[vertex].erase(it);
            break;
        }
    }
}

void BoyerMyrvold::add_to_merge_points(int vertex, bool storeOldHandlesPolicy) {
    mergeStack.push_back(make_tuple(vertex, storeOldHandlesPolicy, false));
}

vector<int> BoyerMyrvold::findExternalFace(int root) {
    vector<int> externalFace;
    set<int> visited;
    stack<int> stack;
    stack.push(root);

    while (!stack.empty()) {
        int v = stack.top();
        stack.pop();

        if (visited.find(v) == visited.end()) {
            visited.insert(v);
            externalFace.push_back(v);

            for (int w : g.adj[v]) {
                if (visited.find(w) == visited.end()) {
                    stack.push(w);
                }
            }
        }
    }

    return externalFace;
}

void BoyerMyrvold::performWalkup(int v) {
    typedef face_vertex_iterator<first_side> walkup_iterator_t;

    // Process all adjacents to find backedges
    for (int w : g.adj[v]) {
        if (w == v) continue; // Skip self-loops

        // Only process backedges to descendants (higher DFS number)
        if (dfsNum[w] <= dfsNum[v] || w == dfsParent[v]) continue;

        // Record backedge
        backedges[v].push_back(make_pair(v, w));
        backedges[w].push_back(make_pair(w, v)); // Also record the other direction

        int timestamp = dfsNum[v];
        backedgeFlag[w] = timestamp; // Mark w as having a backedge to v

        // Initialize walkup iterator starting at w
        walkup_iterator_t walkup_itr(w, &faceHandles);
        walkup_iterator_t walkup_end;

        int lead_vertex = w;

        // Begin walkup traversal
        while (true) {
            // Follow the face until we find a visited vertex
            while (walkup_itr != walkup_end && visited[*walkup_itr] != timestamp) {
                lead_vertex = *walkup_itr;
                visited[lead_vertex] = timestamp;
                ++walkup_itr;
            }

            // Handle end of path
            if (walkup_itr == walkup_end || walkup_itr.currentVertex == -1) {
                int dfs_child = canonicalDfsChild[lead_vertex];
                if (dfs_child != -1) {
                    int parent = dfsParent[dfs_child];
                    if (parent != -1) {
                        // Mark face endpoints as visited
                        int firstVertex = dfsChildHandles[dfs_child].firstVertex();
                        int secondVertex = dfsChildHandles[dfs_child].secondVertex();

                        if (firstVertex >= 0 && firstVertex < visited.size())
                            visited[firstVertex] = timestamp;
                        if (secondVertex >= 0 && secondVertex < visited.size())
                            visited[secondVertex] = timestamp;

                        // Add to pertinentRoots - CRITICAL PART FOR FIXING THE ISSUE
                        if (lowPoint[dfs_child] < dfsNum[v] ||
                            leastAncestor[dfs_child] < dfsNum[v]) {
                            // Higher priority - lower lowPoint or leastAncestor
                            pertinentRoots[parent].push_back(dfsChildHandles[dfs_child]);
                        } else {
                            pertinentRoots[parent].push_front(dfsChildHandles[dfs_child]);
                        }

                        // Continue walkup from parent if not already visited and not at v
                        if (parent != v &&
                            parent >= 0 && parent < visited.size() &&
                            visited[parent] != timestamp) {
                            walkup_itr = walkup_iterator_t(parent, &faceHandles);
                            lead_vertex = parent;
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            } else {
                break;
            }
        }
    }
}

bool BoyerMyrvold::performWalkdown(int v) {
    int w; // the other endpoint of the edge we're embedding

    while (!pertinentRoots[v].empty()) {
        face_handle_t rootFaceHandle = pertinentRoots[v].front();
        face_handle_t currFaceHandle = rootFaceHandle;
        pertinentRoots[v].pop_front();

        mergeStack.clear();

        while (true) {
            auto firstFaceItr = face_vertex_iterator<first_side>(currFaceHandle.firstVertex(), &faceHandles);
            auto secondFaceItr = face_vertex_iterator<second_side>(currFaceHandle.secondVertex(), &faceHandles);
            auto firstFaceEnd = face_vertex_iterator<first_side>(-1, &faceHandles);
            auto secondFaceEnd = face_vertex_iterator<second_side>(-1, &faceHandles);

            int firstSideVertex = -1;
            int secondSideVertex = -1;
            int firstTail = currFaceHandle.getAnchor();
            int secondTail = firstTail;

            for (; firstFaceItr != firstFaceEnd; ++firstFaceItr) {
                int faceVertex = *firstFaceItr;
                if (pertinent(faceVertex, v) || externally_active(faceVertex, v)) {
                    firstSideVertex = faceVertex;
                    secondSideVertex = faceVertex;
                    break;
                }
                firstTail = faceVertex;
            }

            if (firstSideVertex == -1 || firstSideVertex == currFaceHandle.getAnchor())
                break;

            for (; secondFaceItr != secondFaceEnd; ++secondFaceItr) {
                int faceVertex = *secondFaceItr;
                if (pertinent(faceVertex, v) || externally_active(faceVertex, v)) {
                    secondSideVertex = faceVertex;
                    break;
                }
                secondTail = faceVertex;
            }

            int chosen;
            bool choseFirstUpperPath;
            if (internally_active(firstSideVertex, v)) {
                chosen = firstSideVertex;
                choseFirstUpperPath = true;
            } else if (internally_active(secondSideVertex, v)) {
                chosen = secondSideVertex;
                choseFirstUpperPath = false;
            } else if (pertinent(firstSideVertex, v)) {
                chosen = firstSideVertex;
                choseFirstUpperPath = true;
            } else if (pertinent(secondSideVertex, v)) {
                chosen = secondSideVertex;
                choseFirstUpperPath = false;
            } else {
                for (; *firstFaceItr != secondSideVertex; ++firstFaceItr) {
                    int p = *firstFaceItr;
                    if (pertinent(p, v)) {
                        kuratowskiV = v;
                        kuratowskiX = firstSideVertex;
                        kuratowskiY = secondSideVertex;
                        return false;
                    }
                }

                if (firstSideVertex == secondSideVertex) {
                    if (firstTail != v) {
                        int first = faceHandles[firstTail].firstVertex();
                        int second = faceHandles[firstTail].secondVertex();
                        tie(firstSideVertex, firstTail) = make_tuple(firstTail, (first == firstSideVertex) ? second : first);
                    } else if (secondTail != v) {
                        int first = faceHandles[secondTail].firstVertex();
                        int second = faceHandles[secondTail].secondVertex();
                        tie(secondSideVertex, secondTail) = make_tuple(secondTail, (first == secondSideVertex) ? second : first);
                    } else {
                        break;
                    }
                }

                canonicalDfsChild[firstSideVertex] = canonicalDfsChild[rootFaceHandle.firstVertex()];
                canonicalDfsChild[secondSideVertex] = canonicalDfsChild[rootFaceHandle.secondVertex()];
                rootFaceHandle.setFirstVertex(firstSideVertex);
                rootFaceHandle.setSecondVertex(secondSideVertex);

                if (faceHandles[firstSideVertex].firstVertex() == firstTail)
                    faceHandles[firstSideVertex].setFirstVertex(v);
                else
                    faceHandles[firstSideVertex].setSecondVertex(v);

                if (faceHandles[secondSideVertex].firstVertex() == secondTail)
                    faceHandles[secondSideVertex].setFirstVertex(v);
                else
                    faceHandles[secondSideVertex].setSecondVertex(v);

                break;
            }

            bool choseFirstLowerPath = (choseFirstUpperPath && faceHandles[chosen].firstVertex() == firstTail) ||
                                       (!choseFirstUpperPath && faceHandles[chosen].firstVertex() == secondTail);

            if (backedgeFlag[chosen] == dfsNum[v]) {
                w = chosen;
                backedgeFlag[chosen] = g.V + 1; // Mark as processed
                add_to_merge_points(chosen, true);

                for (const auto& e : backedges[chosen]) {
                    add_to_embedded_edges(e, storeOldHandlesPolicy);

                    if (choseFirstLowerPath)
                        faceHandles[chosen].pushFirst(e, g);
                    else
                        faceHandles[chosen].pushSecond(e, g);
                }
            } else {
                mergeStack.push_back(make_tuple(chosen, choseFirstUpperPath, choseFirstLowerPath));
                currFaceHandle = *pertinentRoots[chosen].begin();
                continue;
            }

            bool bottomPathFollowsFirst;
            bool topPathFollowsFirst;
            bool nextBottomFollowsFirst = choseFirstUpperPath;
            face_handle_t topHandle, bottomHandle;
            int mergePoint = chosen;

            while (!mergeStack.empty()) {
                bottomPathFollowsFirst = nextBottomFollowsFirst;
                tie(mergePoint, nextBottomFollowsFirst, topPathFollowsFirst) = mergeStack.back();
                mergeStack.pop_back();

                topHandle = faceHandles[mergePoint];
                bottomHandle = *pertinentRoots[mergePoint].begin();

                int bottomDfsChild = canonicalDfsChild[pertinentRoots[mergePoint].begin()->firstVertex()];
                remove_vertex_from_separated_dfs_child_list(mergePoint);
                pertinentRoots[mergePoint].pop_front();

                add_to_merge_points(topHandle.getAnchor(), true);

                if (topPathFollowsFirst && bottomPathFollowsFirst) {
                    bottomHandle.flip();
                    topHandle.glueFirstToSecond(bottomHandle);
                } else if (!topPathFollowsFirst && bottomPathFollowsFirst) {
                    flipped[bottomDfsChild] = true;
                    topHandle.glueSecondToFirst(bottomHandle);
                } else if (topPathFollowsFirst && !bottomPathFollowsFirst) {
                    flipped[bottomDfsChild] = true;
                    topHandle.glueFirstToSecond(bottomHandle);
                } else {
                    bottomHandle.flip();
                    topHandle.glueSecondToFirst(bottomHandle);
                }
            }

            canonicalDfsChild[w] = canonicalDfsChild[rootFaceHandle.firstVertex()];
            add_to_merge_points(rootFaceHandle.getAnchor(), true);

            for (const auto& e : backedges[chosen]) {
                if (nextBottomFollowsFirst)
                    rootFaceHandle.pushFirst(e, g);
                else
                    rootFaceHandle.pushSecond(e, g);
            }

            backedges[chosen].clear();
            currFaceHandle = rootFaceHandle;
        }
    }

    return true;
}

void BoyerMyrvold::store_old_face_handles(StoreOldHandlesPolicy& policy) {
    // Use the policy object to store old face handles
    for (int v = 0; v < g.V; v++) {
        policy(faceHandles, v);
    }
}

void BoyerMyrvold::clean_up_embedding(StoreEmbeddingPolicy& policy) {
    // Use the policy object to clean up the embedding
    for (int v = 0; v < g.V; v++) {
        policy(faceHandles, v);
    }
}

void BoyerMyrvold::add_to_embedded_edges(const pair<int, int>& e, StoreOldHandlesPolicy& policy) {
    // Implementation for embedding edges
    int u = e.first;
    int v = e.second;

    // Store old handles using the policy before modifying the embedding
    policy(faceHandles, u);
    policy(faceHandles, v);


    // Example implementation:
    if (faceHandles[u].firstVertex() == -1) {
        // Initialize face handle if needed
        faceHandles[u].setFirstVertex(v);
        faceHandles[u].setSecondVertex(u);
    }

    if (faceHandles[v].firstVertex() == -1) {
        // Initialize face handle if needed
        faceHandles[v].setFirstVertex(u);
        faceHandles[v].setSecondVertex(v);
    }
}