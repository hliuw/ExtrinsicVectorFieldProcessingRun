#pragma once

#include <string>

#include "Geometry.h"

namespace dirichlet::io {

bool ReadPlyMeshWithFaceVectors(const std::string& path, TriangleMesh& mesh, std::string& error_message);

}  // namespace dirichlet::io
