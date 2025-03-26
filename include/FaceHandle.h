#ifndef FACEHANDLE_H
#define FACEHANDLE_H

#include <vector>
#include <utility>
#include <iostream>

#include "face_vertex_iterator.h"
#include "Graph.h"

// Forward declare template class to allow friendship
class face_handle_t {
public:
    face_handle_t() : first_vertex(-1), second_vertex(-1), anchor(-1) {}

    int firstVertex() const;
    int secondVertex() const;
    void setFirstVertex(int v);
    void setSecondVertex(int v);
    int getAnchor() const;
    void pushFirst(const std::pair<int, int>& e, const Graph& g);
    void pushSecond(const std::pair<int, int>& e, const Graph& g);
    void flip();
    void glueFirstToSecond(const face_handle_t& other);
    void glueSecondToFirst(const face_handle_t& other);
    void setAnchor(int a);

    // Debug method to check state
    void print() const {
        std::cout << "Face: " << first_vertex << " -> " << second_vertex
                  << " (anchor: " << anchor << ", edges: " << edges.size() << ")" << std::endl;
    }

private:
    int first_vertex;
    int second_vertex;
    int anchor;
    std::vector<std::pair<int, int>> edges;

    // Declare specific instantiations as friends
    friend class face_vertex_iterator<first_side>;
    friend class face_vertex_iterator<second_side>;
    friend class face_vertex_iterator<both_sides>;
};


#endif // FACEHANDLE_H