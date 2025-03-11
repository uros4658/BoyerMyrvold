#ifndef EMBEDDING_H
#define EMBEDDING_H

#include <vector>

struct HalfEdge;
struct Face;

struct HalfEdge {
    int origin;            // Origin vertex
    HalfEdge* twin;        // Twin half-edge
    HalfEdge* next;        // Next half-edge in the face
    int faceID;            // ID of the face this half-edge belongs to

    HalfEdge(int orig);
};

struct Face {
    std::vector<HalfEdge*> boundary; // Half-edges that form the boundary
    Face* prev;                      // Previous face in the list
    Face* next;                      // Next face in the list

    Face();
    void addEdge(HalfEdge* edge);
};

class Embedding {
public:
    std::vector<HalfEdge*> halfEdges;  // All half-edges
    std::vector<Face*> faces;          // All faces
    Face* faceHead;                    // First face in the list
    int dfsTime;                       // For DFS
    std::vector<int> disc;             // Discovery times
    std::vector<int> low;              // Lowpoint values

    Embedding();
    ~Embedding();

    void addEdge(int u, int v);
    void addFace(Face* f);
    void constructFaces();
    void dfsEmbeddingLowpoint(int u, std::vector<bool>& visited, std::vector<HalfEdge*>& lastEdge, const std::vector<std::vector<HalfEdge*>>& outgoing);
    void buildOutgoingMapLowpoint(int numVertices);
};

#endif