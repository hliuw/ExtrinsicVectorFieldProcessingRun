#include <cmath>
#include <iostream>
#include <string>

#include "MeshIO.h"
#include "VectorDirichletEnergy.h"

namespace {

struct CommandLineOptions {
    std::string input_path;
    std::string output_path;
    dirichlet::interpolation::VertexToFaceInterpolationMethod interpolation_method =
        dirichlet::interpolation::VertexToFaceInterpolationMethod::kNone;
    bool use_unfold = true;
    bool show_help = false;
};

dirichlet::TriangleMesh BuildSampleMesh() {
    using dirichlet::Triangle;
    using dirichlet::TriangleMesh;

    TriangleMesh mesh;
    mesh.vertices = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {1.0, 1.0, 1.0}};
    mesh.vertex_colors = {
        {255, 0, 0},
        {0, 255, 0},
        {0, 0, 255},
        {255, 255, 255}};
    mesh.faces = {
        Triangle{{0, 1, 2}},
        Triangle{{1, 3, 2}}};
    mesh.face_vectors = {
        {0.0, 0.0, 2.0},
        {-1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)}};
    return mesh;
}

void PrintUsage(const char* executable_name) {
    std::cout << "Usage: " << executable_name
              << " [--interpolation average|rotation] [--nounfold] [--output output_mesh.ply] [input_mesh.ply]"
              << std::endl;
    std::cout << "  Default input: .ply with face vf_0/vf_1/vf_2 properties" << std::endl;
    std::cout << "  --interpolation, -i   Interpolate vertex vf_* to faces. Supported values: average, rotation"
              << std::endl;
    std::cout << "  --nounfold            Compare neighboring face vectors directly without normal transport"
              << std::endl;
    std::cout << "  --output, -o          Write the loaded/interpolated face-based vf_* mesh to a .ply file" << std::endl;
    std::cout << "  --help, -h            Show this help message" << std::endl;
}

bool ParseCommandLine(int argc, char** argv, CommandLineOptions& options, std::string& error_message) {
    for (int argument_index = 1; argument_index < argc; ++argument_index) {
        const std::string argument = argv[argument_index];

        if (argument == "--help" || argument == "-h") {
            options.show_help = true;
            continue;
        }

        if (argument == "--nounfold") {
            options.use_unfold = false;
            continue;
        }

        if (argument == "--interpolation" || argument == "-i") {
            if (argument_index + 1 >= argc) {
                error_message = "Expected an interpolation method after " + argument + ".";
                return false;
            }

            ++argument_index;
            if (!dirichlet::interpolation::TryParseVertexToFaceInterpolationMethod(
                    argv[argument_index], options.interpolation_method)) {
                error_message = "Unsupported interpolation method '" + std::string(argv[argument_index]) + "'.";
                return false;
            }
            continue;
        }

        const std::string interpolation_prefix = "--interpolation=";
        if (argument.rfind(interpolation_prefix, 0) == 0) {
            const std::string method_name = argument.substr(interpolation_prefix.size());
            if (!dirichlet::interpolation::TryParseVertexToFaceInterpolationMethod(
                    method_name, options.interpolation_method)) {
                error_message = "Unsupported interpolation method '" + method_name + "'.";
                return false;
            }
            continue;
        }

        if (argument == "--output" || argument == "-o") {
            if (argument_index + 1 >= argc) {
                error_message = "Expected an output file path after " + argument + ".";
                return false;
            }

            ++argument_index;
            options.output_path = argv[argument_index];
            continue;
        }

        const std::string output_prefix = "--output=";
        if (argument.rfind(output_prefix, 0) == 0) {
            options.output_path = argument.substr(output_prefix.size());
            continue;
        }

        if (!argument.empty() && argument[0] == '-') {
            error_message = "Unknown command line argument '" + argument + "'.";
            return false;
        }

        if (!options.input_path.empty()) {
            error_message = "Please provide only one input mesh path.";
            return false;
        }

        options.input_path = argument;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    CommandLineOptions options;
    std::string command_line_error;
    if (!ParseCommandLine(argc, argv, options, command_line_error)) {
        std::cerr << command_line_error << std::endl;
        PrintUsage(argv[0]);
        return 1;
    }

    if (options.show_help) {
        PrintUsage(argv[0]);
        return 0;
    }

    dirichlet::TriangleMesh mesh;

    if (!options.input_path.empty()) {
        std::string error_message;
        if (!dirichlet::io::ReadMeshWithFaceVectors(
                options.input_path, mesh, error_message, options.interpolation_method)) {
            std::cerr << "Could not read mesh: " << error_message << std::endl;
            return 1;
        }
    } else {
        mesh = BuildSampleMesh();
        std::cout << "No input mesh provided. Running on a built-in tetrahedral sample mesh." << std::endl;
    }

    if (!options.output_path.empty()) {
        std::string error_message;
        if (!dirichlet::io::WritePlyMeshWithFaceVectors(options.output_path, mesh, error_message)) {
            std::cerr << "Could not write mesh: " << error_message << std::endl;
            return 1;
        }
    }

    std::string validation_error;
    if (!dirichlet::ValidateMeshForDirichletEnergy(mesh, validation_error)) {
        std::cerr << "Invalid mesh for energy evaluation: " << validation_error << std::endl;
        return 1;
    }

    const dirichlet::DirichletEnergyResult result =
        dirichlet::ComputeVectorDirichletEnergy(mesh, options.use_unfold);
    std::cout << "Interior edges: " << result.interior_edge_count << std::endl;
    std::cout << "Total weighted energy: " << result.total_weighted_energy << std::endl;
    std::cout << "Average weighted energy: " << result.average_weighted_energy << std::endl;
    std::cout << "Area weighted vector norm: " << result.area_weighted_vector_norm << std::endl;
    if (result.area_weighted_vector_norm > 0.0) {
        std::cout << "Total energy / vector norm: " << result.total_weighted_energy / result.area_weighted_vector_norm << std::endl;
    }
    return 0;
}
