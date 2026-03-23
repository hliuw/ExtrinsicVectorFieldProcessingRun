#include <iostream>
#include <string>

#include "MeshIO.h"
#include "VectorDirichletEnergy.h"

namespace {

dirichlet::TriangleMesh BuildSampleMesh() {
    using dirichlet::Triangle;
    using dirichlet::TriangleMesh;
    // just for testing
    TriangleMesh mesh;
    mesh.vertices = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {1.0, 1.0, 1.0}};
    mesh.faces = {
        Triangle{{0, 1, 2}},
        Triangle{{1, 3, 2}}};
    mesh.face_vectors = {
        {0.0, 0.0, 1.0},
        {-1/sqrt(3), -1/sqrt(3), 1/sqrt(3)}};
    return mesh;
}

}  // namespace

int main(int argc, char** argv) {
    dirichlet::TriangleMesh mesh;

    if (argc > 1) {
        std::string error_message;
        if (!dirichlet::io::ReadPlyMeshWithFaceVectors(argv[1], mesh, error_message)) {
            std::cerr << "Could not read mesh: " << error_message << std::endl;
            return 1;
        }
    } else {
        mesh = BuildSampleMesh();
        std::cout << "No .ply input provided. Running on a built-in tetrahedral sample mesh." << std::endl;
    }

    std::string validation_error;
    if (!dirichlet::ValidateMeshForDirichletEnergy(mesh, validation_error)) {
        std::cerr << "Invalid mesh for energy evaluation: " << validation_error << std::endl;
        return 1;
    }

    const dirichlet::DirichletEnergyResult result = dirichlet::ComputeVectorDirichletEnergy(mesh);
    std::cout << "Interior edges: " << result.interior_edge_count << std::endl;
    std::cout << "Total weighted energy: " << result.total_weighted_energy << std::endl;
    std::cout << "Average weighted energy: " << result.average_weighted_energy << std::endl;
    return 0;
}
