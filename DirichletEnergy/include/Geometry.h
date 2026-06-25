#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <ostream>
#include <utility>
#include <vector>

namespace dirichlet {

using VertexColor = std::array<unsigned char, 3>;

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Triangle {
    std::array<std::size_t, 3> vertex_indices{};
};

struct TriangleMesh {
    std::vector<Vec3> vertices;
    std::vector<VertexColor> vertex_colors;
    std::vector<Triangle> faces;
    std::vector<Vec3> face_vectors;
};

inline Vec3 operator+(const Vec3& lhs, const Vec3& rhs) {
    return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

inline Vec3 operator-(const Vec3& lhs, const Vec3& rhs) {
    return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

inline Vec3 operator*(const Vec3& value, double scale) {
    return {value.x * scale, value.y * scale, value.z * scale};
}

inline Vec3 operator*(double scale, const Vec3& value) {
    return value * scale;
}

inline Vec3 operator/(const Vec3& value, double scale) {
    return {value.x / scale, value.y / scale, value.z / scale};
}

inline double Dot(const Vec3& lhs, const Vec3& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

inline Vec3 Cross(const Vec3& lhs, const Vec3& rhs) {
    return {
        lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.z * rhs.x - lhs.x * rhs.z,
        lhs.x * rhs.y - lhs.y * rhs.x};
}

inline double SquaredNorm(const Vec3& value) {
    return Dot(value, value);
}

inline double Norm(const Vec3& value) {
    return std::sqrt(SquaredNorm(value));
}

inline Vec3 Normalize(const Vec3& value) {
    const double length = Norm(value);
    if (length == 0.0) {
        return {};
    }
    return value / length;
}

inline Vec3 TriangleNormal(const TriangleMesh& mesh, std::size_t face_index) {
    const Triangle& face = mesh.faces.at(face_index);
    const Vec3& p0 = mesh.vertices.at(face.vertex_indices[0]);
    const Vec3& p1 = mesh.vertices.at(face.vertex_indices[1]);
    const Vec3& p2 = mesh.vertices.at(face.vertex_indices[2]);
    return Normalize(Cross(p1 - p0, p2 - p0));
}

// Returns the number of directed edges that appear in more than one face.
// A consistently wound mesh has each directed edge (a→b) in exactly one face;
// its neighbor uses the opposite direction (b→a). A non-zero return value means
// some adjacent face pairs have inconsistent winding.
inline std::size_t CountWindingInconsistencies(const TriangleMesh& mesh) {
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> directed_edge_count;
    for (const Triangle& face : mesh.faces) {
        for (int i = 0; i < 3; ++i) {
            const std::size_t a = face.vertex_indices[i];
            const std::size_t b = face.vertex_indices[(i + 1) % 3];
            directed_edge_count[{a, b}]++;
        }
    }
    std::size_t inconsistent_count = 0;
    for (const auto& [edge, count] : directed_edge_count) {
        if (count > 1) {
            ++inconsistent_count;
        }
    }
    return inconsistent_count;
}

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

}  // namespace dirichlet
