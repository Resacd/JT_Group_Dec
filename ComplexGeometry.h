#pragma once
#ifndef _MESH_GENERATOR_H_
#define _MESH_GENERATOR_H_

#include "Struct.h"
#include <tuple>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <fstream>
#include "Struct.h"
#include <unordered_set>
#include <mpi.h>
#include <omp.h>
#include <map>
#include "GeometryWriter.h"
#include "GeometryReader.h"
#include "mpi_call.h"
#include "Utility.h"

class Geometry {

public:

	Writer Writer;
	Reader Reader;
	Utility Utility;
	MPI_CALL<MPI_Datatype> mpi_caller;
	std::vector<std::vector<std::vector<Point>>> pointsBox;
	std::vector<Triangle>globalTriangles;
	std::vector<Triangle>localTriangles;
	std::string stlFileName;
	BoundingBox boundingbox;
        std::vector<double> inletOutletCoords;
	std::vector<Point> boundaryPoints;
	std::vector<Point> InletboundaryPoints;
	std::vector<Point> OutletboundaryPoints1;
	std::vector<Point> OutletboundaryPoints2;
	std::vector<Point> OutletboundaryPoints3;
	std::vector<Point> OutletboundaryPoints4;
	std::vector<Point> WallboundaryPoints;
	int globalNumtris;
	std::vector<int> GeomSize;
	double spacing;
	double rayLength;
	int myrank, num_procs;
        int nx,ny,nz;
	std::shared_ptr<VelocitySet> velSet_;

	Geometry() = default;
	Geometry(const std::string& FileName, std::shared_ptr<VelocitySet> velSet,double dx, std::vector<double>& inletOutletCoords);
	~Geometry();
	void GetBoundingbox();
	void DistributeTriangles();
	bool isPointInsideBoundingBox(double x, double y, double z, const BoundingBox& box);
	void generatePointsInBoundingBox();
	bool DoesLineBoxIntersect(BoundingBox& box, Line& line, double t0, double t1);
	bool doLinesIntersect(const Line& line1, const Line& line2);
        bool isPointInsideTriangle(double x, double y, double z, const Triangle& triangle);
        int countIntersectinsTriangle(Line& line, std::vector<Triangle>& posintersectedtri);  
	std::vector<Point> findBoundaryPoints();
	void adjustLineCoordinates(Line& line, BoundingBox& boundingbox, double rayOffset, int directionSwitch);
	void labelPoints();
	std::vector<Point> excludeInletOutlet();
        void JustifyInOutWall();
	void getInOutBP();
};

#endif
