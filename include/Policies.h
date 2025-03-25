#ifndef POLICIES_H
#define POLICIES_H

struct StoreOldHandlesPolicy {
    void operator()(std::vector<face_handle_t>& faceHandles, int vertex) {
        // Implementation of the policy to store old face handles
        face_handle_t oldHandle = faceHandles[vertex];
        // Store or process the old handle as needed
    }
};

struct StoreEmbeddingPolicy {
    void operator()(std::vector<face_handle_t>& faceHandles, int vertex) {
        // Implementation of the policy to clean up embedding
        faceHandles[vertex] = face_handle_t();
        // Perform additional cleanup or processing as needed
    }
};

#endif // POLICIES_H