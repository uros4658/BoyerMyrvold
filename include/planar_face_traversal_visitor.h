class planar_face_traversal_visitor {
public:
    planar_face_traversal_visitor() = default;

    void begin_face() {
        current_face.clear();
        std::cout << "Starting traversal of a new face." << std::endl;
    }

    void next_vertex(int vertex) {
        current_face.push_back(vertex);
        std::cout << "Visited vertex: " << vertex << std::endl;
    }

    void end_face() {
        std::cout << "Finished traversal of face with vertices: ";
        for (int vertex : current_face) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    }

private:
    std::vector<int> current_face; // Stores the vertices of the current face
};