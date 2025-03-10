#include "../include/Embedding.h"
using namespace std;

HalfEdge::HalfEdge(int orig) : origin(orig), twin(nullptr), next(nullptr), faceID(-1) {}

Face::Face() : prev(nullptr), next(nullptr) {}

void Face::addEdge(HalfEdge* edge) {
    boundary.push_back(edge);
    edge->faceID = boundary.size() - 1;
}

Embedding::Embedding() : faceHead(nullptr), dfsTime(0) {}

Embedding::~Embedding() {
    for (auto edge : halfEdges) {
        delete edge;
    }

    Face* current = faceHead;
    while (current != nullptr) {
        Face* next = current->next;
        delete current;
        current = next;
    }
}

void Embedding::addEdge(int u, int v) {
    HalfEdge* e1 = new HalfEdge(u);
    HalfEdge* e2 = new HalfEdge(v);

    // Connect the twin edges
    e1->twin = e2;
    e2->twin = e1;

    halfEdges.push_back(e1);
    halfEdges.push_back(e2);
}

void Embedding::addFace(Face* f) {
    faces.push_back(f);

    if (faceHead == nullptr) {
        faceHead = f;
    } else {
        Face* current = faceHead;
        while (current->next != nullptr) {
            current = current->next;
        }
        current->next = f;
        f->prev = current;
    }
}

void Embedding::constructFaces() {
    // Identify and construct faces from the half-edge structure
    vector<bool> visited(halfEdges.size(), false);

    for (size_t i = 0; i < halfEdges.size(); i++) {
        if (!visited[i]) {
            Face* newFace = new Face();

            // Traverse the boundary of the face
            HalfEdge* start = halfEdges[i];
            HalfEdge* current = start;

            do {
                visited[current - halfEdges[0]] = true;
                newFace->addEdge(current);

                // Find the next edge in the face
                HalfEdge* nextEdge = nullptr;
                for (auto edge : halfEdges) {
                    if (edge->origin == current->twin->origin &&
                        edge != current->twin) {
                        nextEdge = edge;
                        break;
                    }
                }

                if (nextEdge == nullptr) break;
                current = nextEdge;

            } while (current != start && !visited[current - halfEdges[0]]);

            addFace(newFace);
        }
    }
}

void Embedding::buildOutgoingMapLowpoint(int numVertices) {
    vector<vector<HalfEdge*>> outgoing(numVertices);

    // Create outgoing edges map
    for (auto edge : halfEdges) {
        outgoing[edge->origin].push_back(edge);
    }

    // Initialize discovery and lowpoint arrays
    disc.resize(numVertices, -1);
    low.resize(numVertices, -1);
    vector<bool> visited(numVertices, false);
    vector<HalfEdge*> lastEdge(numVertices, nullptr);

    dfsTime = 0;

    // Run DFS to compute lowpoints
    for (int i = 0; i < numVertices; i++) {
        if (!visited[i]) {
            dfsEmbeddingLowpoint(i, visited, lastEdge, outgoing);
        }
    }
}

void Embedding::dfsEmbeddingLowpoint(int u, vector<bool>& visited, vector<HalfEdge*>& lastEdge,
                                   const vector<vector<HalfEdge*>>& outgoing) {
    visited[u] = true;
    disc[u] = low[u] = ++dfsTime;

    for (auto edge : outgoing[u]) {
        int v = edge->twin->origin;

        if (!visited[v]) {
            lastEdge[v] = edge;

            dfsEmbeddingLowpoint(v, visited, lastEdge, outgoing);

            low[u] = min(low[u], low[v]);
        }
        else if (lastEdge[u] != edge->twin) {
            low[u] = min(low[u], disc[v]);
        }
    }
}