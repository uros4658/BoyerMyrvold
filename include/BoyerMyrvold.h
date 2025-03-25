#ifndef BOYERMYRVOLD_H
#define BOYERMYRVOLD_H

#include <list>
#include <vector>
#include <stack>
#include <map>
#include <set>
#include "Graph.h"
#include "FaceHandle.h"
#include "Policies.h"


class BoyerMyrvold {
public:
    BoyerMyrvold(Graph& graph);

    void initializeDataStructures();
    bool isPlanar();
    void performWalkup(int v);
    bool performWalkdown(int v);
    bool pertinent(int vertex, int v);
    bool externally_active(int vertex, int v);
    bool internally_active(int vertex, int v);
    void remove_vertex_from_separated_dfs_child_list(int vertex);
    void add_to_merge_points(int vertex, bool storeOldHandlesPolicy);
    std::vector<int> findExternalFace(int root);

private:
    Graph& g;
    std::vector<bool> visited;
    std::vector<int> dfsNum;
    std::vector<int> backedgeFlag;
    std::vector<int> canonicalDfsChild;
    std::vector<int> dfsParent;
    std::vector<face_handle_t> dfsChildHandles;
    std::vector<int> lowPoint;
    std::vector<int> leastAncestor;
    std::vector<std::deque<face_handle_t>> pertinentRoots;
    std::vector<face_handle_t> faceHandles;
    std::vector<bool> flipped;
    std::vector<int> vertices_by_dfs_num;
    std::vector<std::vector<std::pair<int, int>>> backedges;
    std::stack<std::tuple<int, bool, bool>> mergeStack;
    int kuratowskiV;
    int kuratowskiX;
    int kuratowskiY;
    std::vector<std::pair<int, int>> dfsParentEdge;

    void store_old_face_handles(StoreOldHandlesPolicy& policy);
    void clean_up_embedding(StoreEmbeddingPolicy& policy);
    void add_to_embedded_edges(const std::pair<int, int>& e, StoreOldHandlesPolicy& policy);
};

#endif // BOYERMYRVOLD_H