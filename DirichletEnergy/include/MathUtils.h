#pragma once

#include <algorithm>

#include "Geometry.h"

namespace dirichlet {

inline constexpr double kGeometricEpsilon = 1e-12;

inline double Clamp(double value, double min_value, double max_value) {
    return std::max(min_value, std::min(max_value, value));
}

inline Vec3 AnyOrthogonalUnitVector(const Vec3& normal) {
    Vec3 axis = Cross(normal, {1.0, 0.0, 0.0});
    if (Norm(axis) < kGeometricEpsilon) {
        axis = Cross(normal, {0.0, 1.0, 0.0});
    }
    return Normalize(axis);
}

inline Vec3 RotateAroundAxis(const Vec3& value, const Vec3& axis, double cosine_angle, double sine_angle) {
    return value * cosine_angle + Cross(axis, value) * sine_angle + axis * Dot(axis, value) * (1.0 - cosine_angle);
}

inline Vec3 RotateVectorBetweenNormals(const Vec3& value, const Vec3& from_normal, const Vec3& to_normal) {
    const Vec3 from = Normalize(from_normal);
    const Vec3 to = Normalize(to_normal);
    const double cosine_angle = Clamp(Dot(from, to), -1.0, 1.0);

    if (cosine_angle > 1.0 - kGeometricEpsilon) {
        return value;
    }

    if (cosine_angle < -1.0 + kGeometricEpsilon) {
        const Vec3 axis = AnyOrthogonalUnitVector(from);
        return RotateAroundAxis(value, axis, -1.0, 0.0);
    }

    const Vec3 cross_value = Cross(from, to);
    const Vec3 axis = Normalize(cross_value);
    const double sine_angle = Norm(cross_value);
    return RotateAroundAxis(value, axis, cosine_angle, sine_angle);
}

}  // namespace dirichlet
