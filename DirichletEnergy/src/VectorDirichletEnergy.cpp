#include "VectorDirichletEnergy.h"

#include <algorithm>
#include <array>
#include <functional>
#include <unordered_map>

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

double Clamp(double value, double min_value, double max_value) {
    return std::max(min_value, std::min(max_value, value));
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

Vec3 AnyOrthogonalUnitVector(const Vec3& normal) {
    Vec3 axis = Cross(normal, {1.0, 0.0, 0.0});
    if (Norm(axis) < 1e-12) {
        axis = Cross(normal, {0.0, 1.0, 0.0});
    }
    return Normalize(axis);
}

Vec3 RotateAroundAxis(const Vec3& value, const Vec3& axis, double cosine_angle, double sine_angle) { // Rodrigues' rotation formula
    return value * cosine_angle + Cross(axis, value) * sine_angle + axis * Dot(axis, value) * (1.0 - cosine_angle); 
}

Vec3 RotateVectorBetweenNormals(const Vec3& value, const Vec3& from_normal, const Vec3& to_normal) {
    const Vec3 from = Normalize(from_normal);
    const Vec3 to = Normalize(to_normal);
    const double cosine_angle = Clamp(Dot(from, to), -1.0, 1.0);

    if (cosine_angle > 1.0 - 1e-12) {
        return value;
    }

    if (cosine_angle < -1.0 + 1e-12) {
        const Vec3 axis = AnyOrthogonalUnitVector(from);
        return RotateAroundAxis(value, axis, -1.0, 0.0);
    }

    const Vec3 cross_value = Cross(from, to);
    const Vec3 axis = Normalize(cross_value);
    const double sine_angle = Norm(cross_value);
    return RotateAroundAxis(value, axis, cosine_angle, sine_angle);
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

DirichletEnergyResult ComputeVectorDirichletEnergy(const TriangleMesh& mesh) {
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

        const Vec3 n1 = TriangleNormal(mesh, face1);
        const Vec3 n2 = TriangleNormal(mesh, face2);
        const Vec3 transported_v2 = RotateVectorBetweenNormals(mesh.face_vectors[face2], n2, n1);
        const double mismatch = Norm(mesh.face_vectors[face1] - transported_v2);

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
