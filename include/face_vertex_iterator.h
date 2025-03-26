// face_vertex_iterator.h
#pragma once
#include <vector>

class face_handle_t;

enum Side { first_side, second_side, both_sides };

template <Side side>
class face_vertex_iterator {
public:
    face_vertex_iterator(int startVertex, const std::vector<face_handle_t>& faceHandles)
        : currentVertex(startVertex), faceHandles(faceHandles) {
        if (side == both_sides) {
            // Initialize for both sides
            firstSideItr = faceHandles[startVertex].firstVertex();
            secondSideItr = faceHandles[startVertex].secondVertex();
        } else if (side == first_side) {
            // Initialize for first side
            firstSideItr = faceHandles[startVertex].firstVertex();
        } else {
            // Initialize for second side
            secondSideItr = faceHandles[startVertex].secondVertex();
        }
    }

    int operator*() const {
        if (side == both_sides) {
            return (firstSideItr != -1) ? firstSideItr : secondSideItr;
        } else if (side == first_side) {
            return firstSideItr;
        } else {
            return secondSideItr;
        }
    }

    face_vertex_iterator& operator++() {
        if (side == both_sides) {
            if (firstSideItr != -1) {
                firstSideItr = faceHandles[firstSideItr].firstVertex();
            } else {
                secondSideItr = faceHandles[secondSideItr].secondVertex();
            }
        } else if (side == first_side) {
            firstSideItr = faceHandles[firstSideItr].firstVertex();
        } else {
            secondSideItr = faceHandles[secondSideItr].secondVertex();
        }
        return *this;
    }

    bool operator!=(const face_vertex_iterator& other) const {
        return currentVertex != other.currentVertex;
    }

    int currentVertex;

private:
    const std::vector<face_handle_t>& faceHandles;
    int firstSideItr = -1;
    int secondSideItr = -1;
};