#pragma once

#include <array>
#include <cmath>
#include <ostream>
#include <cstddef>
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

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

}  // namespace dirichlet
