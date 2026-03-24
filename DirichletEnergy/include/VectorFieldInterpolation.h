#pragma once

#include <string>
#include <vector>

#include "Geometry.h"

namespace dirichlet::interpolation {

enum class VertexToFaceInterpolationMethod {
    kNone,
    kVertexAverage,
};

bool TryParseVertexToFaceInterpolationMethod(const std::string& method_name,
                                             VertexToFaceInterpolationMethod& method);
bool InterpolateVertexVectorsToFaceVectors(const TriangleMesh& mesh,
                                           const std::vector<Vec3>& vertex_vectors,
                                           VertexToFaceInterpolationMethod method,
                                           std::vector<Vec3>& face_vectors,
                                           std::string& error_message);

}  // namespace dirichlet::interpolation
