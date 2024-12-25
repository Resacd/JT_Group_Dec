#include "ComplexGeometry.h"
#include <queue>
#include <unordered_set>
#include <map>


Geometry::Geometry(const std::string& FileName, std::shared_ptr<VelocitySet> velSet, double dx,std::vector<double> &inletOutletCoords) {

    stlFileName = FileName;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    velSet_=velSet;
    this->spacing = dx;
    GetBoundingbox();
    this->inletOutletCoords=inletOutletCoords;
    DistributeTriangles();
    generatePointsInBoundingBox();
    labelPoints();
   // Writer.writeLabeledPointsToFile(pointsBox, "LabeledPointsFromGeometry.vtk");
    boundaryPoints = findBoundaryPoints();
   // Writer.writePointsToFile(boundaryPoints, "boundaryPointsFromGeometry.vtk");
    getInOutBP();
   // Writer.writePointsToFile(InletboundaryPoints, "InboundaryPointsFromGeometry.vtk");
   // Writer.writePointsToFile(OutletboundaryPoints1,"OutboundaryPointsFromGeometry1.vtk");
   // Writer.writePointsToFile(OutletboundaryPoints2,"OutboundaryPointsFromGeometry2.vtk");
   // Writer.writePointsToFile(OutletboundaryPoints3,"OutboundaryPointsFromGeometry3.vtk");
   // Writer.writePointsToFile(OutletboundaryPoints4,"OutboundaryPointsFromGeometry4.vtk");
    WallboundaryPoints = excludeInletOutlet();
   // Writer.writePointsToFile(WallboundaryPoints, "WallboundaryPointsFromGeometry.vtk");
}


Geometry::~Geometry() {


}

void Geometry::GetBoundingbox() {


    globalTriangles = Reader.readSTLASCII(stlFileName);
    if (myrank == 0) {

        globalNumtris = globalTriangles.size();
        Point minPoint = { std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
        Point maxPoint = { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };


        for (int t = 0; t < globalTriangles.size(); ++t) {
            const auto& triangle = globalTriangles[t];
            {
                minPoint.x = std::min(minPoint.x, std::min(triangle.x1, std::min(triangle.x2, triangle.x3)));
                minPoint.y = std::min(minPoint.y, std::min(triangle.y1, std::min(triangle.y2, triangle.y3)));
                minPoint.z = std::min(minPoint.z, std::min(triangle.z1, std::min(triangle.z2, triangle.z3)));
                maxPoint.x = std::max(maxPoint.x, std::max(triangle.x1, std::max(triangle.x2, triangle.x3)));
                maxPoint.y = std::max(maxPoint.y, std::max(triangle.y1, std::max(triangle.y2, triangle.y3)));
                maxPoint.z = std::max(maxPoint.z, std::max(triangle.z1, std::max(triangle.z2, triangle.z3)));
            }
        }

        boundingbox = { minPoint.x,maxPoint.x, minPoint.y,maxPoint.y, minPoint.z, maxPoint.z };
    }
    MPI_Bcast(&globalNumtris, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Datatype mpiBoundingBoxType = mpi_caller.createBoundingBoxMPIType();
    MPI_Bcast(&boundingbox, 1, mpiBoundingBoxType, 0, MPI_COMM_WORLD);


    int nx = std::ceil((boundingbox.max_x - boundingbox.min_x) / spacing);
    int ny = std::ceil((boundingbox.max_y - boundingbox.min_y) / spacing);
    int nz = std::ceil((boundingbox.max_z - boundingbox.min_z) / spacing);

    GeomSize = { nx + 1, ny + 1, nz + 1 };
}



void Geometry::generatePointsInBoundingBox() {

     nx = std::ceil((boundingbox.max_x - boundingbox.min_x) / spacing);
     ny = std::ceil((boundingbox.max_y - boundingbox.min_y) / spacing);
     nz = std::ceil((boundingbox.max_z - boundingbox.min_z) / spacing);

    GeomSize = { nx + 1, ny + 1, nz + 1 };

    pointsBox.resize(nx + 1);
    for (int i = 0; i <= nx; ++i) {
        pointsBox[i].resize(ny + 1);
        for (int j = 0; j <= ny; ++j) {
            pointsBox[i][j].resize(nz + 1);
        }
    }

    double x, y, z;

    for (int k = 0; k <= nz; ++k) {
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                x = boundingbox.min_x + i * spacing;
                y = boundingbox.min_y + j * spacing;
                z = boundingbox.min_z + k * spacing;
                Point p{ x, y, z };
                pointsBox[i][j][k] = p;
            }
        }
    }
}


void Geometry::DistributeTriangles() {
    int localNumTriangles = globalNumtris / num_procs;
    int remainingTriangles = globalNumtris % num_procs;
    spacing*=2;
    std::vector<int> sendcounts(num_procs, localNumTriangles);

    for (int i = 0; i < remainingTriangles; ++i) {
        sendcounts[i]++;
    }

    std::vector<int> displacements(num_procs, 0);

    for (int i = 1; i < num_procs; ++i) {
        displacements[i] = displacements[i - 1] + sendcounts[i - 1];
    }


    localTriangles.resize(sendcounts[myrank]);
    MPI_Datatype mpiTriangleType = mpi_caller.createTriangleMPIType();
    MPI_Scatterv(globalTriangles.data(), sendcounts.data(), displacements.data(), mpiTriangleType,
        localTriangles.data(), sendcounts[myrank], mpiTriangleType, 0, MPI_COMM_WORLD);

}



bool Geometry::doLinesIntersect(const Line& line1, const Line& line2) {

    double x1 = line1.start.x;
    double y1 = line1.start.y;
    double z1 = line1.start.z;
    double x2 = line1.end.x;
    double y2 = line1.end.y;
    double z2 = line1.end.z;

    double x3 = line2.start.x;
    double y3 = line2.start.y;
    double z3 = line2.start.z;
    double x4 = line2.end.x;
    double y4 = line2.end.y;
    double z4 = line2.end.z;

    double ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1));
    double ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1));

    const double epsilon = 1e-7;
    if (ua >= -epsilon && ua <= 1 + epsilon && ub >= -epsilon && ub <= 1 + epsilon) {
        return true;
    }

    return false;
}


bool Geometry::isPointInsideTriangle(double x, double y, double z, const Triangle& triangle) {

    double x1 = triangle.x1, y1 = triangle.y1, z1 = triangle.z1;
    double x2 = triangle.x2, y2 = triangle.y2, z2 = triangle.z2;
    double x3 = triangle.x3, y3 = triangle.y3, z3 = triangle.z3;


    double dirX = 0.0, dirY = 0.0, dirZ = 1.0;


    double uX = x2 - x1, uY = y2 - y1, uZ = z2 - z1;
    double vX = x3 - x1, vY = y3 - y1, vZ = z3 - z1;


    double nX = uY * vZ - uZ * vY;
    double nY = uZ * vX - uX * vZ;
    double nZ = uX * vY - uY * vX;


    double nLengthSq = nX * nX + nY * nY + nZ * nZ;
    const double tol = 1e-16;
    if (nLengthSq < tol * tol) return false;


    double w0X = x - x1, w0Y = y - y1, w0Z = z - z1;


    double a = -(nX * w0X + nY * w0Y + nZ * w0Z);
    double b = nX * dirX + nY * dirY + nZ * dirZ;
    if (fabs(b) < tol) return false;


    double pos = a / b;


    double ix = x + pos * dirX;
    double iy = y + pos * dirY;
    double iz = z + pos * dirZ;


    double uu = uX * uX + uY * uY + uZ * uZ;
    double uv = uX * vX + uY * vY + uZ * vZ;
    double vv = vX * vX + vY * vY + vZ * vZ;
    double wX = ix - x1, wY = iy - y1, wZ = iz - z1;
    double wu = wX * uX + wY * uY + wZ * uZ;
    double wv = wX * vX + wY * vY + wZ * vZ;
    double D = uv * uv - uu * vv;

    double s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0) return false;
    double t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0) return false;

    return true;
}


int Geometry::countIntersectinsTriangle(Line& line, std::vector<Triangle>& posintersectedtri) {
    const double threshold = 1e-16;
    int count = 0;

    for (const auto& triangle : posintersectedtri) {
        double px = line.start.x;
        double py = line.start.y;
        double pz = line.start.z;
        double dx = line.end.x - line.start.x;
        double dy = line.end.y - line.start.y;
        double dz = line.end.z - line.start.z;


        double nx = triangle.nx;
        double ny = triangle.ny;
        double nz = triangle.nz;
        double qx = triangle.x1;
        double qy = triangle.y1;
        double qz = triangle.z1;


        double numerator = (qx - px) * nx + (qy - py) * ny + (qz - pz) * nz;
        double denominator = dx * nx + dy * ny + dz * nz;


        if (fabs(denominator) < threshold) {
            continue;
        }

        double t = numerator / denominator;


        if (t >= 0.0 && t <= 1.0) {

            double ix = px + t * dx;
            double iy = py + t * dy;
            double iz = pz + t * dz;

            if (isPointInsideTriangle(ix, iy, iz, triangle)) {
                count++;
            }
        }
    }

    return count;
}



void Geometry::adjustLineCoordinates(Line& line, BoundingBox& boundingbox, double rayOffset, int directionSwitch) {
    switch (directionSwitch) {
    case 0:  // positive z-direction
        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        if (line.start.y == boundingbox.min_y) {
            line.start.y += rayOffset;
            line.end.y += rayOffset;
        }

        if (line.start.y == boundingbox.max_y) {
            line.start.y -= rayOffset;
            line.end.y -= rayOffset;
        }

        break;

    case 1:  // negative z-direction
        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        if (line.start.y == boundingbox.min_y) {
            line.start.y += rayOffset;
            line.end.y += rayOffset;
        }

        if (line.start.y == boundingbox.max_y) {
            line.start.y -= rayOffset;
            line.end.y -= rayOffset;
        }

        break;

    case 2:  // positive x-direction
    case 3:  // negative x-direction
        if (line.start.z == boundingbox.min_z) {
            line.start.z += rayOffset;
            line.end.z += rayOffset;
        }

        if (line.start.z == boundingbox.max_z) {
            line.start.z -= rayOffset;
            line.end.z -= rayOffset;
        }

        if (line.start.y == boundingbox.min_y) {
            line.start.y += rayOffset;
            line.end.y += rayOffset;
        }

        if (line.start.y == boundingbox.max_y) {
            line.start.y -= rayOffset;
            line.end.y -= rayOffset;
        }
        break;

    case 4:  // positive y-direction
    case 5:  // negative y-direction
        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        if (line.start.z == boundingbox.min_z) {
            line.start.z += rayOffset;
            line.end.z += rayOffset;
        }

        if (line.start.z == boundingbox.max_z) {
            line.start.z -= rayOffset;
            line.end.z -= rayOffset;
        }
        break;
    }
}


void Geometry::labelPoints() {
    const double rayOffset = 0.01;
    std::vector<Triangle> possiblelocalIntersectedTriangles;
    int directionSwitch;
    int local_intersected_triangles;
    int total_intersected_triangles;
    int intersectingCount;
    rayLength = boundingbox.max_z - boundingbox.min_z;
    for (int i = 0; i < pointsBox.size(); i++) {
        for (int j = 0; j < pointsBox[i].size(); j++) {
            for (int k = 0; k < pointsBox[i][j].size(); k++) {
                Line line;
                line.start.x = pointsBox[i][j][k].x;
                line.start.y = pointsBox[i][j][k].y;
                line.start.z = pointsBox[i][j][k].z;

                directionSwitch = 0;
                local_intersected_triangles = 0;
                total_intersected_triangles = 0;
                bool foundNonZeroIntersections = false;

  //              do {
                switch (directionSwitch) {
                case 0:  // positive z-direction
                    rayLength = boundingbox.max_z - boundingbox.min_z;
                    line.end.x = pointsBox[i][j][k].x;
                    line.end.y = pointsBox[i][j][k].y;
                    line.end.z = pointsBox[i][j][k].z + rayLength + 0.1;
                    break;
                case 1:  // negative z-direction
                    rayLength = boundingbox.max_z - boundingbox.min_z;
                    line.end.x = pointsBox[i][j][k].x;
                    line.end.y = pointsBox[i][j][k].y;
                    line.end.z = pointsBox[i][j][k].z - (rayLength + 0.1);
                    break;
                case 2:  // positive x-direction
                    rayLength = boundingbox.max_x - boundingbox.min_x;
                    line.end.x = pointsBox[i][j][k].x + rayLength + 0.1;
                    line.end.y = pointsBox[i][j][k].y;
                    line.end.z = pointsBox[i][j][k].z;
                    break;
                case 3:  // negative x-direction
                    rayLength = boundingbox.max_x - boundingbox.min_x;
                    line.end.x = pointsBox[i][j][k].x - (rayLength + 0.1);
                    line.end.y = pointsBox[i][j][k].y;
                    line.end.z = pointsBox[i][j][k].z;
                    break;
                case 4:  // positive y-direction
                    rayLength = boundingbox.max_y - boundingbox.min_y;
                    line.end.x = pointsBox[i][j][k].x;
                    line.end.y = pointsBox[i][j][k].y + rayLength + 0.1;
                    line.end.z = pointsBox[i][j][k].z;
                    break;
                case 5:  // negative y-direction
                    rayLength = boundingbox.max_y - boundingbox.min_y;
                    line.end.x = pointsBox[i][j][k].x;
                    line.end.y = pointsBox[i][j][k].y - (rayLength + 0.1);
                    line.end.z = pointsBox[i][j][k].z;
                    break;
                }

                adjustLineCoordinates(line, boundingbox, rayOffset, directionSwitch);
                local_intersected_triangles = countIntersectinsTriangle(line, localTriangles);
                MPI_Allreduce(&local_intersected_triangles, &total_intersected_triangles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

                if (total_intersected_triangles % 2 == 1) {
                    pointsBox[i][j][k].Type = NodeType::Fluid;
                }
                else {
                    pointsBox[i][j][k].Type = NodeType::NotFluid;
                }
            }
        }
    }
}


std::vector<Point> Geometry::findBoundaryPoints() {
    for (int i = 0; i < pointsBox.size(); i++) {
        for (int j = 0; j < pointsBox[i].size(); j++) {
            for (int k = 0; k < pointsBox[i][j].size(); k++) {
                Point& point = pointsBox[i][j][k];

                if (point.Type == NodeType::Fluid) {
                    bool isBoundaryPoint = false;
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                if (i + dx >= 0 && i + dx < pointsBox.size() &&
                                    j + dy >= 0 && j + dy < pointsBox[i].size() &&
                                    k + dz >= 0 && k + dz < pointsBox[i][j].size()) {
                                    Point& neighbor = pointsBox[i + dx][j + dy][k + dz];
                                    if (neighbor.Type == NodeType::NotFluid) {
                                        isBoundaryPoint = true;
                                        break;
                                    }
                                }
                            }
                            if (isBoundaryPoint) {
                                break;
                            }
                        }
                        if (isBoundaryPoint) {
                            break;
                        }
                    }
                    bool isAtBoundary = (i == 0 || i == pointsBox.size() - 1) ||
                        (j == 0 || j == pointsBox[i].size() - 1) ||
                        (k == 0 || k == pointsBox[i][j].size() - 1);

                    if (isAtBoundary && !isBoundaryPoint && point.Type == NodeType::Fluid) {
                        isBoundaryPoint = true;
                    }


                    if (isBoundaryPoint) {
                        boundaryPoints.push_back(point);
                    }
                }
            }
        }
    }

    return boundaryPoints;
}


std::vector<Point> Geometry::excludeInletOutlet() {
    std::vector<Point> remainingPoints;

    for (const Point& p : boundaryPoints) {
        if (std::find(InletboundaryPoints.begin(), InletboundaryPoints.end(), p) == InletboundaryPoints.end() &&
            std::find(OutletboundaryPoints1.begin(), OutletboundaryPoints1.end(), p) == OutletboundaryPoints1.end() &&
            std::find(OutletboundaryPoints2.begin(), OutletboundaryPoints2.end(), p) == OutletboundaryPoints2.end() &&
            std::find(OutletboundaryPoints3.begin(), OutletboundaryPoints3.end(), p) == OutletboundaryPoints3.end() &&
            std::find(OutletboundaryPoints4.begin(), OutletboundaryPoints4.end(), p) == OutletboundaryPoints4.end()) {
            remainingPoints.push_back(p);
        }
    }

    return remainingPoints;
}



void Geometry::getInOutBP() {
    long double tolerance = 1e-2;

    for (int i = 0; i < boundaryPoints.size(); i++) {
        long double x = boundaryPoints[i].x;
        long double y = boundaryPoints[i].y;
        long double z = boundaryPoints[i].z;

        long double zb_inlet = -172.611L;
        if (std::abs(z - zb_inlet) <= tolerance) {
            InletboundaryPoints.push_back(boundaryPoints[i]);
        }


        long double xb_outlet1 = 14.9714L;
        long double xb_outlet2 = 1.0714L;
        if (std::abs(x - xb_outlet1) <= tolerance) {
            OutletboundaryPoints1.push_back(boundaryPoints[i]);
        }
        else if ((std::abs(x - xb_outlet2) <= tolerance)) {
            OutletboundaryPoints2.push_back(boundaryPoints[i]);
        }

        long double yb_outlet3 = -188.496L;

        if ((std::abs(y - yb_outlet3) <= tolerance) &&
            (1.5 < x && x < 3.1) && (-161.5 < z && z < -160)) {

            OutletboundaryPoints3.push_back(boundaryPoints[i]);
        }

        long double yb_outlet4 = -190.396L;
        if ((std::abs(y - yb_outlet4) <= tolerance) &&
            (6.8 < x && x < 8) && (-157 < z && z < -155)) {
            OutletboundaryPoints4.push_back(boundaryPoints[i]);
        }
    }
}
