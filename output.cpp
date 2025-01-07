#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm> // For std::reverse if needed

// Helper to write binary doubles in big-endian format
void writeBigEndian(std::ofstream& file, double value) {
    char* bytes = reinterpret_cast<char*>(&value);
    std::reverse(bytes, bytes + sizeof(double));
    file.write(bytes, sizeof(double));
}

void writeBinaryVTK(const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Cannot open file for writing!" << std::endl;
        return;
    }

    // Write ASCII header
    file << "# vtk DataFile Version 3.0\n";
    file << "Navier-Stokes Simulation Output\n";
    file << "BINARY\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS 2 2 1\n";
    file << "ORIGIN 0.0 0.0 0.0\n";
    file << "SPACING 1.0 1.0 1.0\n";
    file << "POINT_DATA 4\n";

    // Scalars: u_component
    file << "SCALARS u_component double\n";
    file << "LOOKUP_TABLE default\n";
    std::vector<double> u_component = {1.0, 2.0, 3.0, 4.0};
    for (double value : u_component) {
        writeBigEndian(file, value);
    }

    // Scalars: v_component
    file << "SCALARS v_component double\n";
    file << "LOOKUP_TABLE default\n";
    std::vector<double> v_component = {0.5, 1.5, 2.5, 3.5};
    for (double value : v_component) {
        writeBigEndian(file, value);
    }

    // Scalars: w_component
    file << "SCALARS w_component double\n";
    file << "LOOKUP_TABLE default\n";
    std::vector<double> w_component = {0.0, 0.0, 0.0, 0.0};
    for (double value : w_component) {
        writeBigEndian(file, value);
    }

    // Scalars: pressure
    file << "SCALARS pressure double\n";
    file << "LOOKUP_TABLE default\n";
    std::vector<double> pressure = {1.0, 1.1, 1.2, 1.3};
    for (double value : pressure) {
        writeBigEndian(file, value);
    }

    file.close();
    std::cout << "Binary VTK file written to " << filename << std::endl;
}

int main() {
    writeBinaryVTK("output_binary.vtk");
    return 0;
}
