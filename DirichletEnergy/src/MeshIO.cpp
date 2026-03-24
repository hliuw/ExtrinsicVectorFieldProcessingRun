#include "MeshIO.h"

#include <cstdint>
#include <exception>
#include <string>
#include <utility>
#include <vector>

#include "happly.h"

namespace dirichlet::io {
namespace {

Vec3 ToVec3(const std::array<double, 3>& value) {
    return {value[0], value[1], value[2]};
}

VertexColor ToVertexColor(const std::array<unsigned char, 3>& value) {
    return {value[0], value[1], value[2]};
}

std::array<double, 3> ToArray(const Vec3& value) {
    return {value.x, value.y, value.z};
}

std::array<unsigned char, 3> ToArray(const VertexColor& value) {
    return value;
}

std::vector<VertexColor> DefaultVertexColors(std::size_t count) {
    return std::vector<VertexColor>(count, VertexColor{255, 255, 255});
}

bool LoadVertices(const std::vector<std::array<double, 3>>& vertex_positions,
                  TriangleMesh& mesh,
                  std::string& error_message) {
    if (vertex_positions.empty()) {
        error_message = "PLY file does not contain any vertices.";
        return false;
    }

    mesh.vertices.clear();
    mesh.vertices.reserve(vertex_positions.size());
    for (const std::array<double, 3>& position : vertex_positions) {
        mesh.vertices.push_back(ToVec3(position));
    }
    return true;
}

bool LoadVertexColors(happly::Element& vertex_element,
                      std::size_t expected_vertex_count,
                      TriangleMesh& mesh,
                      std::string& error_message) {
    const bool has_red = vertex_element.hasProperty("red");
    const bool has_green = vertex_element.hasProperty("green");
    const bool has_blue = vertex_element.hasProperty("blue");

    if (!has_red && !has_green && !has_blue) {
        mesh.vertex_colors = DefaultVertexColors(expected_vertex_count);
        return true;
    }

    if (!(has_red && has_green && has_blue)) {
        error_message = "PLY vertex element must either define all of red/green/blue or none of them.";
        return false;
    }

    const std::vector<unsigned char> red = vertex_element.getProperty<unsigned char>("red");
    const std::vector<unsigned char> green = vertex_element.getProperty<unsigned char>("green");
    const std::vector<unsigned char> blue = vertex_element.getProperty<unsigned char>("blue");
    if (red.size() != expected_vertex_count || green.size() != expected_vertex_count || blue.size() != expected_vertex_count) {
        error_message = "PLY vertex color properties do not match the number of vertices.";
        return false;
    }

    mesh.vertex_colors.clear();
    mesh.vertex_colors.reserve(expected_vertex_count);
    for (std::size_t index = 0; index < expected_vertex_count; ++index) {
        mesh.vertex_colors.push_back({red[index], green[index], blue[index]});
    }
    return true;
}

bool LoadFaces(const std::vector<std::vector<std::size_t>>& face_indices,
               TriangleMesh& mesh,
               std::string& error_message) {
    if (face_indices.empty()) {
        error_message = "PLY file does not contain any faces.";
        return false;
    }

    mesh.faces.clear();
    mesh.faces.reserve(face_indices.size());
    for (std::size_t face_index = 0; face_index < face_indices.size(); ++face_index) {
        const std::vector<std::size_t>& indices = face_indices[face_index];
        if (indices.size() != 3) {
            error_message = "Expected every face to have exactly 3 vertex indices, but face " +
                            std::to_string(face_index) + " has " + std::to_string(indices.size()) + ".";
            return false;
        }

        mesh.faces.push_back({{indices[0], indices[1], indices[2]}});
    }
    return true;
}

bool ReadVectorFieldFromElement(happly::Element& element,
                                const std::string& element_name,
                                std::size_t expected_count,
                                std::vector<Vec3>& vector_field,
                                std::string& error_message) {
    for (const std::string& property_name : {std::string("vf_0"), std::string("vf_1"), std::string("vf_2")}) {
        if (!element.hasProperty(property_name)) {
            error_message = "PLY " + element_name + " element is missing required property '" + property_name + "'.";
            return false;
        }
    }

    const std::vector<double> vf_0 = element.getProperty<double>("vf_0");
    const std::vector<double> vf_1 = element.getProperty<double>("vf_1");
    const std::vector<double> vf_2 = element.getProperty<double>("vf_2");

    if (vf_0.size() != expected_count || vf_1.size() != expected_count || vf_2.size() != expected_count) {
        error_message = "PLY " + element_name + " vector properties do not match the expected element count.";
        return false;
    }

    vector_field.clear();
    vector_field.reserve(expected_count);
    for (std::size_t index = 0; index < expected_count; ++index) {
        vector_field.push_back({vf_0[index], vf_1[index], vf_2[index]});
    }
    return true;
}

bool ValidateMeshForWriting(const TriangleMesh& mesh, std::string& error_message) {
    if (mesh.vertices.empty()) {
        error_message = "Cannot write a mesh without vertices.";
        return false;
    }
    if (mesh.faces.empty()) {
        error_message = "Cannot write a mesh without faces.";
        return false;
    }
    if (mesh.face_vectors.size() != mesh.faces.size()) {
        error_message = "Cannot write mesh because the number of face vectors does not match the number of faces.";
        return false;
    }
    return true;
}

}  // namespace

bool ReadMeshWithFaceVectors(const std::string& path,
                             TriangleMesh& mesh,
                             std::string& error_message,
                             interpolation::VertexToFaceInterpolationMethod interpolation_method) {
    return ReadPlyMeshWithFaceVectors(path, mesh, error_message, interpolation_method);
}

bool ReadPlyMeshWithFaceVectors(const std::string& path,
                                TriangleMesh& mesh,
                                std::string& error_message,
                                interpolation::VertexToFaceInterpolationMethod interpolation_method) {
    mesh = {};
    error_message.clear();

    try {
        happly::PLYData ply_data(path);
        happly::Element& vertex_element = ply_data.getElement("vertex");

        const std::vector<std::array<double, 3>> vertex_positions = ply_data.getVertexPositions();
        const std::vector<std::vector<std::size_t>> face_indices = ply_data.getFaceIndices<std::size_t>();

        if (!LoadVertices(vertex_positions, mesh, error_message)) {
            mesh = {};
            return false;
        }
        if (!LoadVertexColors(vertex_element, mesh.vertices.size(), mesh, error_message)) {
            mesh = {};
            return false;
        }
        if (!LoadFaces(face_indices, mesh, error_message)) {
            mesh = {};
            return false;
        }

        if (interpolation_method == interpolation::VertexToFaceInterpolationMethod::kNone) {
            std::vector<Vec3> face_vectors;
            happly::Element& face_element = ply_data.getElement("face");
            if (!ReadVectorFieldFromElement(face_element, "face", mesh.faces.size(), face_vectors, error_message)) {
                mesh = {};
                return false;
            }
            mesh.face_vectors = std::move(face_vectors);
            return true;
        }

        std::vector<Vec3> vertex_vectors;
        if (!ReadVectorFieldFromElement(vertex_element, "vertex", mesh.vertices.size(), vertex_vectors, error_message)) {
            mesh = {};
            return false;
        }
        if (!interpolation::InterpolateVertexVectorsToFaceVectors(
                mesh, vertex_vectors, interpolation_method, mesh.face_vectors, error_message)) {
            mesh = {};
            return false;
        }

        return true;
    } catch (const std::exception& e) {
        mesh = {};
        error_message = e.what();
        return false;
    }
}

bool WritePlyMeshWithFaceVectors(const std::string& path,
                                 const TriangleMesh& mesh,
                                 std::string& error_message) {
    error_message.clear();
    if (!ValidateMeshForWriting(mesh, error_message)) {
        return false;
    }

    try {
        happly::PLYData ply_data;

        std::vector<std::array<double, 3>> vertex_positions;
        vertex_positions.reserve(mesh.vertices.size());
        for (const Vec3& vertex : mesh.vertices) {
            vertex_positions.push_back(ToArray(vertex));
        }
        ply_data.addVertexPositions(vertex_positions);

        std::vector<std::array<unsigned char, 3>> vertex_colors;
        vertex_colors.reserve(mesh.vertices.size());
        if (mesh.vertex_colors.size() == mesh.vertices.size()) {
            for (const VertexColor& color : mesh.vertex_colors) {
                vertex_colors.push_back(ToArray(color));
            }
        } else {
            const std::vector<VertexColor> default_colors = DefaultVertexColors(mesh.vertices.size());
            for (const VertexColor& color : default_colors) {
                vertex_colors.push_back(ToArray(color));
            }
        }
        ply_data.addVertexColors(vertex_colors);

        std::vector<std::vector<uint32_t>> face_indices;
        face_indices.reserve(mesh.faces.size());
        for (const Triangle& face : mesh.faces) {
            face_indices.push_back({
                static_cast<uint32_t>(face.vertex_indices[0]),
                static_cast<uint32_t>(face.vertex_indices[1]),
                static_cast<uint32_t>(face.vertex_indices[2])});
        }
        ply_data.addFaceIndices(face_indices);

        std::vector<double> vf_0;
        std::vector<double> vf_1;
        std::vector<double> vf_2;
        vf_0.reserve(mesh.face_vectors.size());
        vf_1.reserve(mesh.face_vectors.size());
        vf_2.reserve(mesh.face_vectors.size());
        for (const Vec3& vector : mesh.face_vectors) {
            vf_0.push_back(vector.x);
            vf_1.push_back(vector.y);
            vf_2.push_back(vector.z);
        }

        happly::Element& face_element = ply_data.getElement("face");
        face_element.addProperty<double>("vf_0", vf_0);
        face_element.addProperty<double>("vf_1", vf_1);
        face_element.addProperty<double>("vf_2", vf_2);

        ply_data.write(path, happly::DataFormat::Binary);
        return true;
    } catch (const std::exception& e) {
        error_message = e.what();
        return false;
    }
}

}  // namespace dirichlet::io
