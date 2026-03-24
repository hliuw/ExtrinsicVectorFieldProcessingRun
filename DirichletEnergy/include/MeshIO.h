#pragma once

#include <string>

#include "Geometry.h"
#include "VectorFieldInterpolation.h"

namespace dirichlet::io {

bool ReadMeshWithFaceVectors(
    const std::string& path,
    TriangleMesh& mesh,
    std::string& error_message,
    interpolation::VertexToFaceInterpolationMethod interpolation_method =
        interpolation::VertexToFaceInterpolationMethod::kNone);
bool ReadPlyMeshWithFaceVectors(
    const std::string& path,
    TriangleMesh& mesh,
    std::string& error_message,
    interpolation::VertexToFaceInterpolationMethod interpolation_method =
        interpolation::VertexToFaceInterpolationMethod::kNone);
bool WritePlyMeshWithFaceVectors(const std::string& path,
                                 const TriangleMesh& mesh,
                                 std::string& error_message);

}  // namespace dirichlet::io
