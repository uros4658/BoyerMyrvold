#include "../include/FaceHandle.h"
#include "../include/Graph.h"

using namespace std;
int face_handle_t::firstVertex() const {
    return first_vertex;
}

int face_handle_t::secondVertex() const {
    return second_vertex;
}

void face_handle_t::setFirstVertex(int v) {
    first_vertex = v;
}

void face_handle_t::setSecondVertex(int v) {
    second_vertex = v;
}

int face_handle_t::getAnchor() const {
    return anchor;
}

void face_handle_t::pushFirst(const pair<int, int>& e, const Graph& g) {
    edges.push_back(e);
}

void face_handle_t::pushSecond(const pair<int, int>& e, const Graph& g) {
    edges.push_back(e);
}

void face_handle_t::flip() {
    swap(first_vertex, second_vertex);
}

void face_handle_t::glueFirstToSecond(const face_handle_t& other) {
    // Glue the first side of this face to the second side of the other face
    // This combines two faces by connecting their boundaries

    // Store current values
    int old_first = this->first_vertex;
    int old_second = other.second_vertex;

    // Update vertices to maintain face connectivity
    this->setFirstVertex(old_second);

    // Merge the edge lists, ensuring proper orientation
    vector<pair<int, int>> new_edges = this->edges;

    // Add edges from other face, preserving orientation
    for (const auto& edge : other.edges) {
        new_edges.push_back(edge);
    }

    this->edges = new_edges;

    // Update anchor if needed
    if (this->anchor == -1) {
        this->anchor = other.anchor;
    }
}

void face_handle_t::glueSecondToFirst(const face_handle_t& other) {
    // Glue the second side of this face to the first side of the other face
    // This combines two faces by connecting their boundaries

    // Store current values
    int old_second = this->second_vertex;
    int old_first = other.first_vertex;

    // Update vertices to maintain face connectivity
    this->setSecondVertex(old_first);

    // Merge the edge lists, ensuring proper orientation
    vector<pair<int, int>> new_edges = this->edges;

    // Add edges from other face, preserving orientation
    for (const auto& edge : other.edges) {
        new_edges.push_back(edge);
    }

    this->edges = new_edges;

    // Update anchor if needed
    if (this->anchor == -1) {
        this->anchor = other.anchor;
    }
}

void face_handle_t::setAnchor(int a) {
    anchor = a;
}