#include "GeometryWriter.h"


Writer::Writer() {


}

Writer::~Writer() {


}

void Writer::writePointoffToFile(const std::vector<Point>& points, const std::string& filename) {
    // Open the output file for writing
    std::ofstream outFile(filename);

    // Write the header information for the VTK file format
    outFile << "# vtk DataFile Version 3.0" << std::endl;
    outFile << "Points data" << std::endl;
    outFile << "ASCII" << std::endl;
    outFile << "DATASET POLYDATA" << std::endl;

    // Write the points data
    outFile << "POINTS " << points.size() << " double" << std::endl;
    for (const auto& point : points) {
        outFile << point.x << " " << point.y << " " << point.z << std::endl;
    }

    // Close the output file
    outFile.close();
}


void Writer::writePointsToFile(const std::vector<Point>& points, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Write the header information
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "Points Inside Complex Geometry" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write the POINTS section
    int numPoints = points.size();
    file << "POINTS " << numPoints << " double" << std::endl;
    for (const auto& point : points) {
        file << point.x << " " << point.y << " " << point.z << std::endl;
    }

    // Write the CELLS section
    file << "CELLS " << numPoints << " " << numPoints * 2 << std::endl;
    for (int i = 0; i < numPoints; ++i) {
        file << "1 " << i << std::endl;
    }

    // Write the CELL_TYPES section
    file << "CELL_TYPES " << numPoints << std::endl;
    for (int i = 0; i < numPoints; ++i) {
        file << "1" << std::endl;
    }

    file.close();
}



void Writer::writeBoxToFile(const std::vector<std::vector<std::vector<Point>>>& points, const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Unable to open file " << filename << " for writing" << std::endl;
        return;
    }

    int nx = points.size();
    int ny = points[0].size();
    int nz = points[0][0].size();

    int numPoints = nx * ny * nz;

    file << "# vtk DataFile Version 3.0\n";
    file << "Points\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << numPoints << " float\n";

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                file << (points[i][j][k].x) << " " << (points[i][j][k].y) << " " << (points[i][j][k].z) << "\n";
            }
        }
    }


    file.close();
    std::cout << " " << "Points written to file " << filename << std::endl;
}

void Writer::writePointsToFile(const std::vector<std::vector<std::vector<Point>>>& pointsBox, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Write the header information
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "Points Inside Complex Geometry" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Collect the fluid points
    std::vector<Point> fluidPoints;
    for (const auto& plane : pointsBox) {
        for (const auto& row : plane) {
            for (const auto& point : row) {

                    fluidPoints.push_back(point);
            }
        }
    }

    // Write the POINTS section
    file << "POINTS " << fluidPoints.size() << " double" << std::endl;
    for (const auto& point : fluidPoints) {
        file << point.x << " " << point.y << " " << point.z << std::endl;
    }

    // Write the CELLS section
    int numCells = fluidPoints.size();
    file << "CELLS " << numCells << " " << numCells * 2 << std::endl;
    for (int i = 0; i < numCells; ++i) {
        file << "1 " << i << std::endl;
    }

    // Write the CELL_TYPES section
    file << "CELL_TYPES " << numCells << std::endl;
    for (int i = 0; i < numCells; ++i) {
        file << "1" << std::endl;
    }

    file.close();
}


void Writer::writeLabeledPointsToFile(const std::vector<std::vector<std::vector<Point>>>& pointsBox, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Write the header information
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "Points Inside Complex Geometry" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Collect the fluid points
    std::vector<Point> fluidPoints;
    for (const auto& plane : pointsBox) {
        for (const auto& row : plane) {
            for (const auto& point : row) {
                if (point.Type == NodeType::Fluid) {
                    fluidPoints.push_back(point);
                }
            }
        }
    }

    // Write the POINTS section
    file << "POINTS " << fluidPoints.size() << " double" << std::endl;
    for (const auto& point : fluidPoints) {
        file << point.x << " " << point.y << " " << point.z << std::endl;
    }

    // Write the CELLS section
    int numCells = fluidPoints.size();
    file << "CELLS " << numCells << " " << numCells * 2 << std::endl;
    for (int i = 0; i < numCells; ++i) {
        file << "1 " << i << std::endl;
    }

    // Write the CELL_TYPES section
    file << "CELL_TYPES " << numCells << std::endl;
    for (int i = 0; i < numCells; ++i) {
        file << "1" << std::endl;
    }

    file.close();
}

void Writer::writeLineToVTK(const Line& line, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Write the header information
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "Line for Visualization" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write the POINTS section
    file << "POINTS 2 double" << std::endl;
    file << line.start.x << " " << line.start.y << " " << line.start.z << std::endl;
    file << line.end.x << " " << line.end.y << " " << line.end.z << std::endl;

    // Write the CELLS section for the line
    file << "CELLS 1 3" << std::endl;
    file << "2 0 1" << std::endl;

    // Write the CELL_TYPES section (use VTK_LINE)
    file << "CELL_TYPES 1" << std::endl;
    file << "3" << std::endl;

    file.close();
}

void Writer::WriteLatticeToVTK(const std::string& filename, const std::vector<std::shared_ptr<Node>>& lattice,int start_x,
    int start_y, int start_z,double spacing) {
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Unable to open VTK file for writing." << std::endl;
        return;
    }

    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "Lattice Data" << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    vtkFile << "POINTS " << lattice.size() << " double" << std::endl;
    for (const auto& nodePtr : lattice) {
        vtkFile << (nodePtr->x_position + start_x) * spacing << " " << (nodePtr->y_position + start_y) * spacing << " " << (nodePtr->z_position + start_z) * spacing << std::endl;
    }

    vtkFile << "CELLS " << lattice.size() << " " << 2 * lattice.size() << std::endl;
    for (size_t i = 0; i < lattice.size(); ++i) {
        vtkFile << "1 " << i << std::endl;
    }


    vtkFile << "CELL_TYPES " << lattice.size() << std::endl;
    for (size_t i = 0; i < lattice.size(); ++i) {
        vtkFile << "1" << std::endl;
    }

    vtkFile.close();
    std::cout << "VTK file " << filename << "' written successfully." << std::endl;
}


void Writer::WriteNotNullLatToVTKFile(const std::string& filename, std::vector<std::vector<std::vector<std::shared_ptr<Node>>>> &lattice) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return;
    }


    outfile << "# vtk DataFile Version 3.0" << std::endl;
    outfile << "Lattice Nodes" << std::endl;
    outfile << "ASCII" << std::endl;
    outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    int my_Nx_ = lattice.size();
    int my_Ny_ = lattice[0].size();
    int my_Nz_ = lattice[0][0].size();
    int numNonNullNodes = 0;
    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (lattice[i][j][k] != nullptr) {
                    numNonNullNodes++;
                }
            }
        }
    }

    outfile << "POINTS " << numNonNullNodes << " int" << std::endl;
    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (lattice[i][j][k] != nullptr) {
                    outfile << i << " " << j << " " << k << std::endl;
                }
            }
        }
    }

    outfile << "CELLS " << numNonNullNodes << " " << numNonNullNodes * 2 << std::endl;
    for (int cellId = 0; cellId < numNonNullNodes; cellId++) {
        outfile << "1 " << cellId << std::endl;
    }

    outfile << "CELL_TYPES " << numNonNullNodes << std::endl;
    for (int cellId = 0; cellId < numNonNullNodes; cellId++) {
        outfile << "1" << std::endl;
    }

    outfile.close();
}

