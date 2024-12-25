#pragma once
#ifndef _SUBDOMAIN_H_
#define _SUBDOMAIN_H_

#include <vector>
#include "Node.h"
#include "HaloLayer.h"
#include "mpi.h"
#include <omp.h>
#include "typedef.h"
#include "Struct.h"
#include "Boundary.h"
#include "Utility.h"
#include "ComplexGeometry.h"


class SubDomain {

private:
	int myrank_, num_procs_;
	void getNeighborInfo(int neighbor_rank, Neighbors& neighbor);
	bool is_valid_coordinate(int coords[]);
	int  get_rank(MPI_Comm comm, int coords[]);
	
public:

	Writer Writer;
	Utility Utility_;
	Geometry Geometry_;
	int start_x_;
	int	start_y_;
	int	start_z_;
	int	end_x_;
	int	end_y_;
	int	end_z_;
	int my_Nx_;
	int my_Ny_;
	int my_Nz_;
	int halo_size;
	std::vector<int> gridSize_;
	std::vector<int> NewgridSize_;
	std::vector<int> m_subdomainSize;
	int SubdomainSize_;
	std::shared_ptr<VelocitySet> m_velSet;
	std::vector< std::vector< std::vector<std::shared_ptr<Node>>>> m_lattice;
	std::vector<std::shared_ptr<Node>> HaloLattice_;
	std::vector<std::shared_ptr<Node>> Involved_Lat_;
	std::vector<std::shared_ptr<Node>> InletBoundaryNode_;
	std::vector<std::shared_ptr<Node>> OutletBoundaryNode1_;
	std::vector<std::shared_ptr<Node>> OutletBoundaryNode2_;
	std::vector<std::shared_ptr<Node>> OutletBoundaryNode3_;
	std::vector<std::shared_ptr<Node>> OutletBoundaryNode4_;
	std::vector<std::shared_ptr<Node>> WallBoundaryNode_;
	std::vector<std::shared_ptr<Node>> nextToTheWall;
	MPI_Comm cart_comm_;
	int my_up_, my_down_, my_left_, my_right_, my_front_, my_back_;
	int right_up_, right_down_, left_up_, left_down_, up_front_, up_back_;
	int down_front_, down_back_, right_front_, right_back_, left_front_, left_back_;
	std::vector<int>coords_;
	std::vector<int>dims_;
	std::vector<Neighbors> neighbors_;
	std::vector<int>need_communication_;

	std::vector<std::pair<std::tuple<int, int, int>, std::tuple<int, int, int>>> local2global_;
	SubDomain() = default;
	SubDomain(Geometry& Geometry);
	~SubDomain();
	void Create_cart_comm();
	void CreateNodes();
	void Connectivity();
	void IdentifyNextToWallLattices();
	void calculateBoundaryNormals();
	void Initialize();
	bool IsInSubDomain(int coord_x, int coord_y, int coord_z);
	std::vector<int> Local2Global(int coord_x, int coord_y, int coord_z);
	bool IsFluidInSubDomain(int coord_x, int coord_y, int coord_z);
	bool IsHaloLayer(int coordX, int coordY, int coordZ);
	void Need_Comm();
	void GetFluidHaloLat();
	const std::shared_ptr<VelocitySet> GetVelSet();


	friend class Simulation;
	friend class Node;

};
#endif
