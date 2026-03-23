#pragma once

#include <cstddef>
#include <string>

#include "Geometry.h"

namespace dirichlet {

struct DirichletEnergyResult {
    double total_weighted_energy = 0.0;
    double average_weighted_energy = 0.0;
    std::size_t interior_edge_count = 0;
};

bool ValidateMeshForDirichletEnergy(const TriangleMesh& mesh, std::string& error_message);
DirichletEnergyResult ComputeVectorDirichletEnergy(const TriangleMesh& mesh);

}  // namespace dirichlet
