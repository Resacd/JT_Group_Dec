#include "Subdomain.h"

SubDomain::SubDomain(Geometry& Geometry) {

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs_);
    Geometry_ = Geometry;
    m_velSet = Geometry_.velSet_;
    gridSize_ = Geometry_.GeomSize;
    Create_cart_comm();
    CreateNodes();
    Connectivity();
//    IdentifyNextToWallLattices();
    calculateBoundaryNormals();
    Initialize();
    GetFluidHaloLat();
    MPI_Comm_free(&cart_comm_);
    
}


SubDomain::~SubDomain() {

}



void SubDomain::getNeighborInfo(int neighbor_rank, Neighbors& neighbor) {
    neighbor.rank = neighbor_rank;
    if (neighbor_rank != MPI_PROC_NULL && neighbor_rank >= 0) {
        MPI_Cart_coords(cart_comm_, neighbor_rank, 3, neighbor.coords);
    }
    else {
        neighbor.coords[0] = -1; 
        neighbor.coords[1] = -1; 
        neighbor.coords[2] = -1; 
    }
}




bool SubDomain::is_valid_coordinate(int coords[]) {
    return coords[0] >= 0 && coords[0] < dims_[0] && coords[1] >= 0 && coords[1] < dims_[1] && coords[2] >= 0 && coords[2] < dims_[2];
}

int SubDomain::get_rank(MPI_Comm comm, int coords[]) {
    int rank;
    MPI_Cart_rank(comm, coords, &rank);
    return rank;
}

void SubDomain::Create_cart_comm() {
    int const ndims = 3;
    dims_ = { 0, 0, 0 };
    int periods[ndims] = { 0, 0, 0 };
    int reorder = 1;
    coords_ = { 0, 0, 0 };
    int remain_dims[3] = { 1, 1, 1 };
    int Nx = gridSize_[0];
    int Ny = gridSize_[1];
    int Nz = gridSize_[2];
    MPI_Dims_create(num_procs_, ndims, dims_.data());
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims_.data(), periods, reorder, &cart_comm_);

    MPI_Cart_shift(cart_comm_, 0, 1, &my_left_, &my_right_);
    MPI_Cart_shift(cart_comm_, 1, 1, &my_front_, &my_back_);
    MPI_Cart_shift(cart_comm_, 2, 1, &my_down_, &my_up_);
    MPI_Comm_rank(cart_comm_, &myrank_);
    MPI_Cart_coords(cart_comm_, myrank_, ndims, coords_.data());
    int cart_coords[ndims];
    MPI_Cart_get(cart_comm_, ndims, dims_.data(), periods, cart_coords);

    neighbors_.resize(18);

    getNeighborInfo(my_left_, neighbors_[0]);   
    getNeighborInfo(my_right_, neighbors_[1]);  
    getNeighborInfo(my_front_, neighbors_[2]);  
    getNeighborInfo(my_back_, neighbors_[3]);   
    getNeighborInfo(my_up_, neighbors_[4]);     
    getNeighborInfo(my_down_, neighbors_[5]);   

    int right_up_coords[3] = { coords_[0] + 1, coords_[1], coords_[2] + 1 };
    int right_down_coords[3] = { coords_[0] + 1, coords_[1], coords_[2] - 1 };
    int left_up_coords[3] = { coords_[0] - 1, coords_[1], coords_[2] + 1 };
    int left_down_coords[3] = { coords_[0] - 1, coords_[1], coords_[2] - 1 };
    int up_front_coords[3] = { coords_[0], coords_[1] - 1, coords_[2] + 1 };
    int up_back_coords[3] = { coords_[0], coords_[1] + 1, coords_[2] + 1 };
    int down_front_coords[3] = { coords_[0], coords_[1] - 1, coords_[2] - 1 };
    int down_back_coords[3] = { coords_[0], coords_[1] + 1, coords_[2] - 1 };
    int right_front_coords[3] = { coords_[0] + 1, coords_[1] - 1, coords_[2] };
    int right_back_coords[3] = { coords_[0] + 1, coords_[1] + 1, coords_[2] };
    int left_front_coords[3] = { coords_[0] - 1, coords_[1] - 1, coords_[2] };
    int left_back_coords[3] = { coords_[0] - 1, coords_[1] + 1, coords_[2] };

    right_up_ = is_valid_coordinate(right_up_coords) ? get_rank(cart_comm_, right_up_coords) : -1;
    right_down_ = is_valid_coordinate(right_down_coords) ? get_rank(cart_comm_, right_down_coords) : -1;
    left_up_ = is_valid_coordinate(left_up_coords) ? get_rank(cart_comm_, left_up_coords) : -1;
    left_down_ = is_valid_coordinate(left_down_coords) ? get_rank(cart_comm_, left_down_coords) : -1;
    up_front_ = is_valid_coordinate(up_front_coords) ? get_rank(cart_comm_, up_front_coords) : -1;
    up_back_ = is_valid_coordinate(up_back_coords) ? get_rank(cart_comm_, up_back_coords) : -1;
    down_front_ = is_valid_coordinate(down_front_coords) ? get_rank(cart_comm_, down_front_coords) : -1;
    down_back_ = is_valid_coordinate(down_back_coords) ? get_rank(cart_comm_, down_back_coords) : -1;
    right_front_ = is_valid_coordinate(right_front_coords) ? get_rank(cart_comm_, right_front_coords) : -1;
    right_back_ = is_valid_coordinate(right_back_coords) ? get_rank(cart_comm_, right_back_coords) : -1;
    left_front_ = is_valid_coordinate(left_front_coords) ? get_rank(cart_comm_, left_front_coords) : -1;
    left_back_ = is_valid_coordinate(left_back_coords) ? get_rank(cart_comm_, left_back_coords) : -1;


    getNeighborInfo(right_up_, neighbors_[6]);     
    getNeighborInfo(right_down_, neighbors_[7]);   
    getNeighborInfo(left_up_, neighbors_[8]);      
    getNeighborInfo(left_down_, neighbors_[9]);    
    getNeighborInfo(up_front_, neighbors_[10]);    
    getNeighborInfo(up_back_, neighbors_[11]);     
    getNeighborInfo(down_front_, neighbors_[12]);  
    getNeighborInfo(down_back_, neighbors_[13]);   
    getNeighborInfo(right_front_, neighbors_[14]); 
    getNeighborInfo(right_back_, neighbors_[15]);  
    getNeighborInfo(left_front_, neighbors_[16]);  
    getNeighborInfo(left_back_, neighbors_[17]);   

    halo_size = 1;
    my_Nx_ = Nx / dims_[0];
    my_Ny_ = Ny / dims_[1];
    my_Nz_ = Nz / dims_[2];
    start_x_ = coords_[0] * my_Nx_;
    start_y_ = coords_[1] * my_Ny_;
    start_z_ = coords_[2] * my_Nz_;
    end_x_ = start_x_ + my_Nx_ - 1;
    end_y_ = start_y_ + my_Ny_ - 1;
    end_z_ = start_z_ + my_Nz_ - 1;

    if (Nx % dims_[0] > 0 && coords_[0] == dims_[0] - 1) {
        my_Nx_ += Nx % dims_[0];
    }
    if (Ny % dims_[1] > 0 && coords_[1] == dims_[1] - 1) {
        my_Ny_ += Ny % dims_[1];
    }
    if (Nz % dims_[2] > 0 && coords_[2] == dims_[2] - 1) {
        my_Nz_ += Nz % dims_[2];
    }


    if (dims_[0] > 1) {
        if (coords_[0] == 0) {
            my_Nx_ += halo_size;
            end_x_ = start_x_ + my_Nx_ - 1;
        }
        else if (coords_[0] == dims_[0] - 1) {
            my_Nx_ += halo_size;
            start_x_ = start_x_ + 2 * halo_size * coords_[0] - halo_size;
            end_x_ = start_x_ + my_Nx_ - 1;
        }
        else {
            my_Nx_ += 2 * halo_size;
            start_x_ = start_x_ + 2 * halo_size * coords_[0] - halo_size;
            end_x_ = start_x_ + my_Nx_ - 1;
        }
    }

    if (dims_[1] > 1) {
        if (coords_[1] == 0) {
            my_Ny_ += halo_size;
            end_y_ = start_y_ + my_Ny_ - 1;
        }
        else if (coords_[1] == dims_[1] - 1) {
            my_Ny_ += halo_size;
            start_y_ = start_y_ + 2 * halo_size * coords_[1] - halo_size;
            end_y_ = start_y_ + my_Ny_ - 1;
        }
        else {
            my_Ny_ += 2 * halo_size;
            start_y_ = start_y_ + 2 * halo_size * coords_[1] - halo_size;
            end_y_ = start_y_ + my_Ny_ - 1;
        }
    }

    if (dims_[2] > 1) {
        if (coords_[2] == 0) {
            my_Nz_ += halo_size;
            end_z_ = start_z_ + my_Nz_ - 1;
        }
        else if (coords_[2] == dims_[2] - 1) {
            my_Nz_ += halo_size;
            start_z_ = start_z_ + 2 * halo_size * coords_[2] - halo_size;
            end_z_ = start_z_ + my_Nz_ - 1;
        }
        else {
            my_Nz_ += 2 * halo_size;
            start_z_ = start_z_ + 2 * halo_size * coords_[2] - halo_size;
            end_z_ = start_z_ + my_Nz_ - 1;
        }
    }

    NewgridSize_.resize(3);

    NewgridSize_[0] = (dims_[0] - 1) * 2 * halo_size + gridSize_[0];
    NewgridSize_[1] = (dims_[1] - 1) * 2 * halo_size + gridSize_[1];
    NewgridSize_[2] = (dims_[2] - 1) * 2 * halo_size + gridSize_[2];

}


bool SubDomain::IsInSubDomain(int coord_x, int coord_y, int coord_z) {

    if (coord_x >= 0 && coord_x < my_Nx_
        && coord_y >= 0 && coord_y < my_Ny_
        && coord_z >= 0 && coord_z < my_Nz_)
        return true;
    else
        return false;

}

std::vector<int>SubDomain::Local2Global(int coord_x, int coord_y, int coord_z) {

    std::vector<int> L2G;
    int global_coordX = coord_x;
    int global_coordY = coord_y;
    int global_coordZ = coord_z;

    if (dims_[0] > 1) {
        if (coords_[0] == 0) {
            global_coordX = coord_x;
        }
        else if (coords_[0] == dims_[0] - 1) {

            global_coordX = start_x_ + coord_x - 2 * coords_[0] * halo_size;
        }
        else {

            global_coordX = start_x_ + coord_x - 2 * coords_[0] * halo_size;

        }
    }


    if (dims_[1] > 1) {
        if (coords_[1] == 0) {
            ;
            global_coordY = coord_y;
        }
        else if (coords_[1] == dims_[1] - 1) {

            global_coordY = start_y_ + coord_y - 2 * coords_[1] * halo_size;
        }
        else {

            global_coordY = start_y_ + coord_y - 2 * coords_[1] * halo_size;

        }
    }

    if (dims_[2] > 1) {
        if (coords_[2] == 0) {
            ;
            global_coordZ = coord_z;
        }
        else if (coords_[2] == dims_[2] - 1) {

            global_coordZ = start_z_ + coord_z - 2 * coords_[2] * halo_size;
        }
        else {

            global_coordZ = start_z_ + coord_z - 2 * coords_[2] * halo_size;

        }
    }

    L2G = { global_coordX,global_coordY,global_coordZ };
    return L2G;

}




bool SubDomain::IsHaloLayer(int coordX, int coordY, int coordZ) {

    if (my_left_ != MPI_PROC_NULL && coordX == 0) {
        return true;
    }

    if (my_right_ != MPI_PROC_NULL && coordX == my_Nx_ - 1) {
        return true;
    }

    if (my_front_ != MPI_PROC_NULL && coordY == 0) {
        return true;
    }

    if (my_back_ != MPI_PROC_NULL && coordY == my_Ny_ - 1) {
        return true;
    }

    if (my_up_ != MPI_PROC_NULL && coordZ == my_Nz_ - 1) {
        return true;
    }

    if (my_down_ != MPI_PROC_NULL && coordZ == 0) {
        return true;
    }

    return false;
}

bool SubDomain::IsFluidInSubDomain(int coord_x, int coord_y, int coord_z) {
    std::vector<int> L2G = Local2Global(coord_x, coord_y, coord_z);

    if (Geometry_.pointsBox[L2G[0]][L2G[1]][L2G[2]].Type == NodeType::Fluid) {

        return true;
    }
    else {
        return false;
    }
}



void SubDomain::CreateNodes() {

    std::vector<Point> InletBoundary = Geometry_.InletboundaryPoints;
    std::vector<Point> OutletBoundary1 = Geometry_.OutletboundaryPoints1;
    std::vector<Point> OutletBoundary2 = Geometry_.OutletboundaryPoints2;
    std::vector<Point> OutletBoundary3 = Geometry_.OutletboundaryPoints3;
    std::vector<Point> OutletBoundary4 = Geometry_.OutletboundaryPoints4;
    std::vector<Point> WallBoundary = Geometry_.WallboundaryPoints;


    m_lattice.resize(my_Nx_);

    for (int i = 0; i < my_Nx_; i++) {
        m_lattice[i].resize(my_Ny_);
        for (int j = 0; j < my_Ny_; j++) {
            m_lattice[i][j].resize(my_Nz_);
            for (int k = 0; k < my_Nz_; k++) {
                m_lattice[i][j][k] = nullptr;
            }
        }
    }

    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (IsHaloLayer(i, j, k)) {
                    m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Halo);
                }
                else if (IsFluidInSubDomain(i, j, k)) {

                    std::vector<int> L2GIndex = Local2Global(i, j, k);
                    int global_coordX = L2GIndex[0];
                    int global_coordY = L2GIndex[1];
                    int global_coordZ = L2GIndex[2];

                    Point& GlobalcurrentPoint = Geometry_.pointsBox[global_coordX][global_coordY][global_coordZ];
                    if (Utility_.isPointInBoundary(GlobalcurrentPoint, InletBoundary)) {

                        m_lattice[i][j][k] = std::make_shared<Boundary>(i, j, k, m_velSet, NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Inlet;
                        InletBoundaryNode_.push_back(m_lattice[i][j][k]);
                    }

                    else if (Utility_.isPointInBoundary(GlobalcurrentPoint, OutletBoundary1)) {

                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Outlet;
                        OutletBoundaryNode1_.push_back(m_lattice[i][j][k]);
                    }

                    else if (Utility_.isPointInBoundary(GlobalcurrentPoint, OutletBoundary2)) {

                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Outlet;
                        OutletBoundaryNode2_.push_back(m_lattice[i][j][k]);
                    }

                    else if (Utility_.isPointInBoundary(GlobalcurrentPoint, OutletBoundary3)) {

                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Outlet;
                        OutletBoundaryNode3_.push_back(m_lattice[i][j][k]);
                    }

                    else if (Utility_.isPointInBoundary(GlobalcurrentPoint, OutletBoundary4)) {

                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Outlet;
                        OutletBoundaryNode4_.push_back(m_lattice[i][j][k]);
                    }


                    else if (Utility_.isPointInBoundary(GlobalcurrentPoint, WallBoundary)) {

                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k,NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Wall;
                        WallBoundaryNode_.push_back(m_lattice[i][j][k]);
                    }

                    else {
                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Fluid);
                    }
                }
            }
        }
    }

// std::cout << "Processor " << myrank_ << " WallBoundaryNode_ size: " << WallBoundaryNode_.size() << std::endl;

}


void SubDomain::Connectivity() {

    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (m_lattice[i][j][k] != nullptr) {

                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        int x_neighbor = i + m_velSet->GetDirection(dir)[0];
                        int y_neighbor = j + m_velSet->GetDirection(dir)[1];
                        int z_neighbor = k + m_velSet->GetDirection(dir)[2];

                        if (IsInSubDomain(x_neighbor, y_neighbor, z_neighbor) &&
                            m_lattice[x_neighbor][y_neighbor][z_neighbor] != nullptr) {
                            m_lattice[i][j][k]->m_neighbours[dir] = m_lattice[x_neighbor][y_neighbor][z_neighbor];
                        }
                        else {
                            m_lattice[i][j][k]->m_neighbours[dir] = nullptr;
                        }
                    }
                }
            }
        }
    }
}



void SubDomain::IdentifyNextToWallLattices() {
    std::unordered_set<std::shared_ptr<Node>, Node::Hash, Node::Equal> uniqueNextToWall;

    for (auto& wallNode : WallBoundaryNode_) {
        for (int dir = 0; dir < _nLatNodes; dir++) {
            auto neighbor = wallNode->m_neighbours[dir];

            if (neighbor) {
                if (neighbor->NodeType_ == NodeType::Fluid && neighbor->BoundaryType_!=BoundaryType::Wall) {

                    uniqueNextToWall.insert(neighbor);
                }
            }
        }
    }

    nextToTheWall.assign(uniqueNextToWall.begin(), uniqueNextToWall.end());

}


/*

void SubDomain::calculateBoundaryNormals() {


    for (auto& wall_lat : nextToTheWall) {

        wall_lat->normal.resize(_dim);
        std::fill(wall_lat->normal.begin(), wall_lat->normal.end(), 0.0);

        int x_pos = wall_lat->x_position;
        int y_pos = wall_lat->y_position;
        int z_pos = wall_lat->z_position;

        std::vector<int> ret = Local2Global(x_pos, y_pos, z_pos);


        double x = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]].x;
        double y = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]].y;
        double z = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]].z;


        for (const auto& triangle : Geometry_.globalTriangles) {

            double ux = triangle.x2 - triangle.x1;
            double uy = triangle.y2 - triangle.y1;
            double uz = triangle.z2 - triangle.z1;

            double vx = triangle.x3 - triangle.x1;
            double vy = triangle.y3 - triangle.y1;
            double vz = triangle.z3 - triangle.z1;

            double nx = uy * vz - uz * vy;
            double ny = uz * vx - ux * vz;
            double nz = ux * vy - uy * vx;

            double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (norm > 0) {
                nx /= norm;
                ny /= norm;
                nz /= norm;
            }

            double area = 0.5 * norm;

            double centroidX = (triangle.x1 + triangle.x2 + triangle.x3) / 3.0;
            double centroidY = (triangle.y1 + triangle.y2 + triangle.y3) / 3.0;
            double centroidZ = (triangle.z1 + triangle.z2 + triangle.z3) / 3.0;

            double distance = std::sqrt(
                std::pow(x - centroidX, 2) +
                std::pow(y - centroidY, 2) +
                std::pow(z - centroidZ, 2)
            );


            double weight = area / (distance * distance);

            wall_lat->normal[0] += nx * weight;
            wall_lat->normal[1] += ny * weight;
            wall_lat->normal[2] += nz * weight;
        }

   
        double norm = std::sqrt(
            wall_lat->normal[0] * wall_lat->normal[0] +
            wall_lat->normal[1] * wall_lat->normal[1] +
            wall_lat->normal[2] * wall_lat->normal[2]
        );

        if (norm > 0) {
            wall_lat->normal[0] /= norm;
            wall_lat->normal[1] /= norm;
            wall_lat->normal[2] /= norm;
        }
    }
}

*/




void SubDomain::calculateBoundaryNormals() {

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < WallBoundaryNode_.size(); ++i) {
        auto& wall_lat = WallBoundaryNode_[i];


        wall_lat->normal.resize(_dim);
        std::fill(wall_lat->normal.begin(), wall_lat->normal.end(), 0.0);

        int x_pos = wall_lat->x_position;
        int y_pos = wall_lat->y_position;
        int z_pos = wall_lat->z_position;

        std::vector<int> ret = Local2Global(x_pos, y_pos, z_pos);

        double x = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]].x;
        double y = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]].y;
        double z = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]].z;

        for (const auto& triangle : Geometry_.globalTriangles) {
            double ux = triangle.x2 - triangle.x1;
            double uy = triangle.y2 - triangle.y1;
            double uz = triangle.z2 - triangle.z1;

            double vx = triangle.x3 - triangle.x1;
            double vy = triangle.y3 - triangle.y1;
            double vz = triangle.z3 - triangle.z1;

            double nx = uy * vz - uz * vy;
            double ny = uz * vx - ux * vz;
            double nz = ux * vy - uy * vx;

            double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (norm > 0) {
                nx /= norm;
                ny /= norm;
                nz /= norm;
            }

            double area = 0.5 * norm;

            double centroidX = (triangle.x1 + triangle.x2 + triangle.x3) / 3.0;
            double centroidY = (triangle.y1 + triangle.y2 + triangle.y3) / 3.0;
            double centroidZ = (triangle.z1 + triangle.z2 + triangle.z3) / 3.0;

            double distance = std::sqrt(
                std::pow(x - centroidX, 2) +
                std::pow(y - centroidY, 2) +
                std::pow(z - centroidZ, 2)
            );

            double weight = area / (distance * distance);


            wall_lat->normal[0] += nx * weight;
            wall_lat->normal[1] += ny * weight;
            wall_lat->normal[2] += nz * weight;
        }


        double norm = std::sqrt(
            wall_lat->normal[0] * wall_lat->normal[0] +
            wall_lat->normal[1] * wall_lat->normal[1] +
            wall_lat->normal[2] * wall_lat->normal[2]
        );

        if (norm > 0) {
            wall_lat->normal[0] /= norm;
            wall_lat->normal[1] /= norm;
            wall_lat->normal[2] /= norm;
        }
    }
}





void SubDomain::Initialize()
{
    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {

                if (m_lattice[i][j][k] != nullptr) {

                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        m_lattice[i][j][k]->Initialize(m_velSet, dir);
                    }

                }
            }
        }
    }
}



void SubDomain::GetFluidHaloLat() {

    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {

                if (m_lattice[i][j][k] != nullptr) {
                    Involved_Lat_.push_back(m_lattice[i][j][k]);

                }
            }
        }
    }
}

const std::shared_ptr<VelocitySet> SubDomain::GetVelSet()
{
    return m_velSet;
}

