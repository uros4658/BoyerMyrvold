#ifndef FACEHANDLE_H
#define FACEHANDLE_H

#include <vector>
#include <utility>
#include <iostream>
#include "Graph.h"

// Forward declare template class to allow friendship
template<int Side> class face_vertex_iterator;

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
    friend class face_vertex_iterator<0>;
    friend class face_vertex_iterator<1>;
    friend class face_vertex_iterator<2>;
};

const int both_sides = 2;
const int first_side = 0;
const int second_side = 1;

template<int Side>
class face_vertex_iterator {
public:
    face_vertex_iterator(int vertex, const std::vector<face_handle_t>& faceHandles)
        : currentVertex(vertex), faceHandles(faceHandles) {}

    face_vertex_iterator& operator++();
    int operator*() const;
    bool operator!=(const face_vertex_iterator& other) const;
    face_vertex_iterator& operator=(const face_vertex_iterator& other);

    int currentVertex;

private:
    const std::vector<face_handle_t>& faceHandles;
};

template<int Side>
face_vertex_iterator<Side>& face_vertex_iterator<Side>::operator++() {
    if (currentVertex == -1) {
        return *this; // Already at end
    }

    bool found = false;

    if constexpr (Side == first_side) {
        // Look for a face handle where current vertex is the first vertex
        for (size_t i = 0; i < faceHandles.size(); i++) {
            if (faceHandles[i].first_vertex == currentVertex) {
                currentVertex = faceHandles[i].second_vertex;
                found = true;
                break;
            }
        }
    }
    else if constexpr (Side == second_side) {
        // Look for a face handle where current vertex is the second vertex
        for (size_t i = 0; i < faceHandles.size(); i++) {
            if (faceHandles[i].second_vertex == currentVertex) {
                currentVertex = faceHandles[i].first_vertex;
                found = true;
                break;
            }
        }
    }
    else { // both_sides
        // Look for a face handle where current vertex is either first or second vertex
        for (size_t i = 0; i < faceHandles.size(); i++) {
            if (faceHandles[i].first_vertex == currentVertex) {
                currentVertex = faceHandles[i].second_vertex;
                found = true;
                break;
            }
            else if (faceHandles[i].second_vertex == currentVertex) {
                currentVertex = faceHandles[i].first_vertex;
                found = true;
                break;
            }
        }
    }

    if (!found) {
        currentVertex = -1; // End of iteration
    }

    return *this;
}

template<int Side>
int face_vertex_iterator<Side>::operator*() const {
    return currentVertex;
}

template<int Side>
bool face_vertex_iterator<Side>::operator!=(const face_vertex_iterator& other) const {
    return currentVertex != other.currentVertex;
}

template<int Side>
face_vertex_iterator<Side>& face_vertex_iterator<Side>::operator=(const face_vertex_iterator& other) {
    if (this != &other) {
        currentVertex = other.currentVertex;
    }
    return *this;
}

#endif // FACEHANDLE_H