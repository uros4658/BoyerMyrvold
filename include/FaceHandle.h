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

    face_vertex_iterator& operator++() {
        if (currentVertex == -1) {
            return *this; // Already at end
        }

        bool found = false;

        if constexpr (Side == first_side) {
            for (const auto& handle : faceHandles) {
                if (handle.first_vertex == currentVertex) {
                    currentVertex = handle.second_vertex;
                    found = true;
                    break;
                }
            }
        }
        else if constexpr (Side == second_side) {
            for (const auto& handle : faceHandles) {
                if (handle.second_vertex == currentVertex) {
                    currentVertex = handle.first_vertex;
                    found = true;
                    break;
                }
            }
        }
        else { // both_sides
            for (const auto& handle : faceHandles) {
                if (handle.first_vertex == currentVertex) {
                    currentVertex = handle.second_vertex;
                    found = true;
                    break;
                }
                else if (handle.second_vertex == currentVertex) {
                    currentVertex = handle.first_vertex;
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

    int operator*() const {
        return currentVertex;
    }

    bool operator!=(const face_vertex_iterator& other) const {
        return currentVertex != other.currentVertex;
    }

    face_vertex_iterator& operator=(const face_vertex_iterator& other) {
        if (this != &other) {
            currentVertex = other.currentVertex;
        }
        return *this;
    }

    int currentVertex;

private:
    const std::vector<face_handle_t>& faceHandles;
};

#endif // FACEHANDLE_H