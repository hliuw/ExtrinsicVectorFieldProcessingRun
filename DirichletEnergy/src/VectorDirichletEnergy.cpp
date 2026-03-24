#include "VectorDirichletEnergy.h"

#include <array>
#include <functional>
#include <unordered_map>

#include "MathUtils.h"

namespace dirichlet {
namespace {

struct EdgeKey {
    std::size_t a = 0;
    std::size_t b = 0;

    bool operator==(const EdgeKey& other) const {
        return a == other.a && b == other.b;
    }
};

struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& key) const {
        return std::hash<std::size_t>{}(key.a) ^ (std::hash<std::size_t>{}(key.b) << 1U);
    }
};

struct EdgeFaceRecord {
    std::size_t face_index = 0;
    std::size_t opposite_vertex = 0;
};

struct EdgeAdjacency {
    std::array<EdgeFaceRecord, 2> incident_faces{};
    std::size_t count = 0;
};

EdgeKey MakeEdgeKey(std::size_t i, std::size_t j) {
    if (i < j) {
        return {i, j};
    }
    return {j, i};
}

double CotangentAtOppositeVertex(const TriangleMesh& mesh,
                                 const EdgeKey& edge,
                                 std::size_t opposite_vertex_index) {
    const Vec3& opposite = mesh.vertices.at(opposite_vertex_index);
    const Vec3& edge_a = mesh.vertices.at(edge.a);
    const Vec3& edge_b = mesh.vertices.at(edge.b);
    const Vec3 u = edge_a - opposite;
    const Vec3 v = edge_b - opposite;
    const double area_factor = Norm(Cross(u, v));
    if (area_factor == 0.0) {
        return 0.0;
    }
    return Dot(u, v) / area_factor;
}

std::unordered_map<EdgeKey, EdgeAdjacency, EdgeKeyHash> BuildEdgeAdjacency(const TriangleMesh& mesh) {
    std::unordered_map<EdgeKey, EdgeAdjacency, EdgeKeyHash> adjacency;

    for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
        const Triangle& face = mesh.faces[face_index];
        const std::array<std::size_t, 3>& v = face.vertex_indices;
        const std::array<EdgeKey, 3> edges = {MakeEdgeKey(v[0], v[1]), MakeEdgeKey(v[1], v[2]), MakeEdgeKey(v[2], v[0])};
        const std::array<std::size_t, 3> opposite = {v[2], v[0], v[1]};

        for (std::size_t i = 0; i < edges.size(); ++i) {
            EdgeAdjacency& entry = adjacency[edges[i]];
            if (entry.count < entry.incident_faces.size()) {
                entry.incident_faces[entry.count++] = {face_index, opposite[i]};
            }
        }
    }

    return adjacency;
}

}  // namespace

bool ValidateMeshForDirichletEnergy(const TriangleMesh& mesh, std::string& error_message) {
    if (mesh.faces.empty()) {
        error_message = "Mesh has no triangles.";
        return false;
    }

    if (mesh.faces.size() != mesh.face_vectors.size()) {
        error_message = "Expected one vector field entry per triangle.";
        return false;
    }

    for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
        const Triangle& face = mesh.faces[face_index];
        for (std::size_t vertex_index : face.vertex_indices) {
            if (vertex_index >= mesh.vertices.size()) {
                error_message = "Triangle references an out-of-range vertex index.";
                return false;
            }
        }

        if (Norm(TriangleNormal(mesh, face_index)) == 0.0) {
            error_message = "Mesh contains a degenerate triangle.";
            return false;
        }
    }

    // disable closed mesh check
    // const auto adjacency = BuildEdgeAdjacency(mesh);
    // for (const auto& pair : adjacency) {
    //     const EdgeAdjacency& entry = pair.second;
    //     if (entry.count != 2) {
    //         error_message = "Expected a closed manifold mesh with exactly two incident triangles per edge.";
    //         return false;
    //     }
    // }

    error_message.clear();
    return true;
}

DirichletEnergyResult ComputeVectorDirichletEnergy(const TriangleMesh& mesh, bool use_unfold) {
    const auto adjacency = BuildEdgeAdjacency(mesh);
    DirichletEnergyResult result;

    for (const auto& pair : adjacency) {
        const EdgeKey& edge = pair.first;
        const EdgeAdjacency& entry = pair.second;
        if (entry.count != 2) {
            continue;
        }

        const std::size_t face1 = entry.incident_faces[0].face_index;
        const std::size_t face2 = entry.incident_faces[1].face_index;
        const double cot_alpha = CotangentAtOppositeVertex(mesh, edge, entry.incident_faces[0].opposite_vertex);
        const double cot_beta = CotangentAtOppositeVertex(mesh, edge, entry.incident_faces[1].opposite_vertex);
        const double weight = cot_alpha + cot_beta;

        const Vec3 face1_vector = mesh.face_vectors[face1];
        const Vec3 face2_vector = mesh.face_vectors[face2];
        Vec3 comparison_vector = face2_vector;
        
        if (use_unfold) {
            const Vec3 n1 = TriangleNormal(mesh, face1);
            const Vec3 n2 = TriangleNormal(mesh, face2);
            // check consistency of normal direction 
            // if(Dot(n1, n2) <= 0.0)
            // {
            //     std::cout << "Warning: Angle larger than 90 degrees detected between faces " << face1 << " " << n1 << " and " << face2 << " " << n2 << "." << std::endl;
            // }
            comparison_vector = RotateVectorBetweenNormals(face2_vector, n2, n1);
        }
        const double mismatch = Norm(face1_vector - comparison_vector);

        result.total_weighted_energy += weight * mismatch;
        ++result.interior_edge_count;
    }

    if (result.interior_edge_count > 0) {
        result.average_weighted_energy =
            result.total_weighted_energy / static_cast<double>(result.interior_edge_count);
    }

    return result;
}

}  // namespace dirichlet
