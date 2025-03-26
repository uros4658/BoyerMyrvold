// face_vertex_iterator.h
#pragma once
#include <vector>
#include <stdexcept>

class face_handle_t;

enum Side { first_side, second_side, both_sides };

template <Side side>
class face_vertex_iterator {
public:
    // Constructor
    face_vertex_iterator(int startVertex, const std::vector<face_handle_t>* faceHandles)
        : currentVertex(startVertex), faceHandles(faceHandles) {
        if (startVertex == -1 || !faceHandles) {
            // End iterator or null
            firstSideItr = -1;
            secondSideItr = -1;
            return;
        }

        // Initialize based on side
        if (side == first_side || side == both_sides) {
            firstSideItr = startVertex >= 0 && startVertex < faceHandles->size() ?
                          (*faceHandles)[startVertex].firstVertex() : -1;
        } else {
            firstSideItr = -1;
        }

        if (side == second_side || side == both_sides) {
            secondSideItr = startVertex >= 0 && startVertex < faceHandles->size() ?
                           (*faceHandles)[startVertex].secondVertex() : -1;
        } else {
            secondSideItr = -1;
        }
    }

    // Copy constructor
    face_vertex_iterator(const face_vertex_iterator& other)
        : currentVertex(other.currentVertex), faceHandles(other.faceHandles),
          firstSideItr(other.firstSideItr), secondSideItr(other.secondSideItr) {}

    // Assignment operator
    face_vertex_iterator& operator=(const face_vertex_iterator& other) {
        if (this != &other) {
            currentVertex = other.currentVertex;
            faceHandles = other.faceHandles;
            firstSideItr = other.firstSideItr;
            secondSideItr = other.secondSideItr;
        }
        return *this;
    }

    // Dereference
    int operator*() const {
        if (currentVertex == -1) {
            throw std::runtime_error("Attempted to dereference end iterator");
        }
        return currentVertex;
    }

    // Pre-increment
    face_vertex_iterator& operator++() {
        if (currentVertex == -1 || !faceHandles) {
            return *this;  // Already at end
        }

        if (side == first_side) {
            if (firstSideItr >= 0 && firstSideItr < faceHandles->size()) {
                firstSideItr = (*faceHandles)[firstSideItr].firstVertex();
                currentVertex = firstSideItr;
            } else {
                firstSideItr = -1;
                currentVertex = -1;
            }
        } else if (side == second_side) {
            if (secondSideItr >= 0 && secondSideItr < faceHandles->size()) {
                secondSideItr = (*faceHandles)[secondSideItr].secondVertex();
                currentVertex = secondSideItr;
            } else {
                secondSideItr = -1;
                currentVertex = -1;
            }
        } else if (side == both_sides) {
            // Try first side
            if (firstSideItr >= 0 && firstSideItr < faceHandles->size()) {
                int next = (*faceHandles)[firstSideItr].firstVertex();
                if (next >= 0 && next < faceHandles->size()) {
                    firstSideItr = next;
                    currentVertex = firstSideItr;
                    return *this;
                } else {
                    firstSideItr = -1;
                }
            }

            // Try second side if first failed
            if (secondSideItr >= 0 && secondSideItr < faceHandles->size()) {
                int next = (*faceHandles)[secondSideItr].secondVertex();
                if (next >= 0 && next < faceHandles->size()) {
                    secondSideItr = next;
                    currentVertex = secondSideItr;
                    return *this;
                } else {
                    secondSideItr = -1;
                }
            }

            // Both sides done
            currentVertex = -1;
        }

        return *this;
    }

    // Comparison operators
    bool operator==(const face_vertex_iterator& other) const {
        return currentVertex == other.currentVertex;
    }

    bool operator!=(const face_vertex_iterator& other) const {
        return !(*this == other);
    }

    // Public member for access
    int currentVertex;

private:
    const std::vector<face_handle_t>* faceHandles; // Changed to pointer from reference
    int firstSideItr = -1;
    int secondSideItr = -1;
};