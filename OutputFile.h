#pragma once
#include "Node.h"
#include <fstream>
#include <string>
#include <mpi.h>
#include "Utility.h"
#include <filesystem>
#include <memory> 
#include <map>
#include <vector>
#include <tuple>
#include <cmath>
#include <filesystem>
#include <sys/stat.h>
#include <limits>




#ifndef _WIN32
#include <sys/stat.h>
#else
#include <windows.h>
#endif



#include "Struct.h"
#include "HaloLayer.h"
#include "Subdomain.h"

#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkDelaunay3D.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkTriangleFilter.h>
#include <vtkHexahedron.h>
#include <vtkCellArray.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkTriangle.h>
#include <vtkPolygon.h>
#include <vtkIdTypeArray.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkCellArray.h>
#include <vtkVertex.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkDelaunay3D.h>
#include <vtkGeometryFilter.h>
#include <vtkPointData.h>
#include <vtkContourFilter.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkPCANormalEstimation.h>
#include <vtkFillHolesFilter.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>




class OutputFile
{
private:

	std::mutex m_fileMutex;
	int m_outputIdx;
	std::mutex m_dataMutex;
	Utility Utility_;

public:

	OutputFile() = default;
	OutputFile(int outputIdx);
	~OutputFile();

        double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2);
        double calculateMeanNeighboringValues(const std::vector<DataPoint>& data, const DataPoint& point, double searchRadius, char attribute);
        void smoothData(std::vector<DataPoint>& data, double neighborSearchRadius);
        void printElapsedTime(std::chrono::steady_clock::time_point startTime);
        void processFlowData(const std::string& inputFileName, const std::string& outputFileName, double neighborSearchRadius);
	void WriteVTK_LocalData(SubDomain &subdomain);
        void WriteVTK_Collective(SubDomain& subdomain);
	void WriteSurfaceAttribute(SubDomain &subdomain,int iter);

};

