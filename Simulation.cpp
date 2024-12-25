#include "Simulation.h"

Simulation::Simulation(SubDomain& SubDomain,double tau,double pressure,double del_t)
	: SubDomain_(SubDomain)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs_);
    SetRelaxation(tau);
    pressure_ = pressure;
    m_dt = del_t;
    m_omega=1.0/m_relaxation;
    InitializeShearStressAccumulators();
    Iterate();

}

Simulation::~Simulation()
{
}


void Simulation::InitializeShearStressAccumulators()
{
    for (auto& node : SubDomain_.WallBoundaryNode_) {
        node->accumulated_shear_stress = { 0.0, 0.0, 0.0 };
        node->accumulated_shear_stress_magnitude = 0.0;
    }
}






void Simulation::Collision()
{
    for (int f = 0; f < SubDomain_.Involved_Lat_.size(); f++) {

        int i = SubDomain_.Involved_Lat_[f]->x_position;
        int j = SubDomain_.Involved_Lat_[f]->y_position;
        int k = SubDomain_.Involved_Lat_[f]->z_position;

        SubDomain_.m_lattice[i][j][k]->density = 0.0;

        for (int dir = 0; dir < _nLatNodes; dir++) {
            SubDomain_.m_lattice[i][j][k]->density += SubDomain_.m_lattice[i][j][k]->m_distributions[dir];
        }

        
        if (SubDomain_.Involved_Lat_[f]->BoundaryType_ == BoundaryType::Wall) {
            SubDomain_.Involved_Lat_[f]->velocity[0] = 0.0;
            SubDomain_.Involved_Lat_[f]->velocity[1] = 0.0;
            SubDomain_.Involved_Lat_[f]->velocity[2] = 0.0;
        } else {

           
#ifdef D3Q19
            SubDomain_.m_lattice[i][j][k]->velocity[0] =
                ((SubDomain_.m_lattice[i][j][k]->m_distributions[0] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[6] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[7] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[8] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[9]) -
                 (SubDomain_.m_lattice[i][j][k]->m_distributions[1] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[10] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[11] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[12] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[13])) / SubDomain_.m_lattice[i][j][k]->density;

            SubDomain_.m_lattice[i][j][k]->velocity[1] =
                ((SubDomain_.m_lattice[i][j][k]->m_distributions[2] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[6] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[10] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[14] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[15]) -
                 (SubDomain_.m_lattice[i][j][k]->m_distributions[3] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[7] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[11] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[16] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[17])) / SubDomain_.m_lattice[i][j][k]->density;

            SubDomain_.m_lattice[i][j][k]->velocity[2] =
                ((SubDomain_.m_lattice[i][j][k]->m_distributions[4] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[8] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[12] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[14] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[16]) -
                 (SubDomain_.m_lattice[i][j][k]->m_distributions[5] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[9] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[13] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[15] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[17])) / SubDomain_.m_lattice[i][j][k]->density;
#endif

           
#ifdef D3Q27
            SubDomain_.m_lattice[i][j][k]->velocity[0] =
                ((SubDomain_.m_lattice[i][j][k]->m_distributions[1] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[7] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[8] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[11] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[13] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[19] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[21] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[23] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[25]) -
                 (SubDomain_.m_lattice[i][j][k]->m_distributions[2] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[9] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[10] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[12] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[14] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[20] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[22] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[24] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[26])) / SubDomain_.m_lattice[i][j][k]->density;

            SubDomain_.m_lattice[i][j][k]->velocity[1] =
                ((SubDomain_.m_lattice[i][j][k]->m_distributions[3] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[7] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[9] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[15] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[16] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[19] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[21] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[24] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[26]) -
                 (SubDomain_.m_lattice[i][j][k]->m_distributions[4] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[8] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[10] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[17] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[18] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[20] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[22] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[23] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[25])) / SubDomain_.m_lattice[i][j][k]->density;

            SubDomain_.m_lattice[i][j][k]->velocity[2] =
                ((SubDomain_.m_lattice[i][j][k]->m_distributions[5] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[11] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[12] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[15] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[17] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[19] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[22] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[23] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[26]) -
                 (SubDomain_.m_lattice[i][j][k]->m_distributions[6] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[13] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[14] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[16] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[18] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[20] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[21] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[24] +
                  SubDomain_.m_lattice[i][j][k]->m_distributions[25])) / SubDomain_.m_lattice[i][j][k]->density;
#endif
        }


        if (SubDomain_.Involved_Lat_[f]->NodeType_ != NodeType::Halo) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                double du = SubDomain_.m_lattice[i][j][k]->velocity[0] * SubDomain_.m_velSet->GetDirection(dir)[0] +
                            SubDomain_.m_lattice[i][j][k]->velocity[1] * SubDomain_.m_velSet->GetDirection(dir)[1] +
                            SubDomain_.m_lattice[i][j][k]->velocity[2] * SubDomain_.m_velSet->GetDirection(dir)[2];

                double u_sqr = SubDomain_.m_lattice[i][j][k]->velocity[0] * SubDomain_.m_lattice[i][j][k]->velocity[0] +
                               SubDomain_.m_lattice[i][j][k]->velocity[1] * SubDomain_.m_lattice[i][j][k]->velocity[1] +
                               SubDomain_.m_lattice[i][j][k]->velocity[2] * SubDomain_.m_lattice[i][j][k]->velocity[2];

                SubDomain_.m_lattice[i][j][k]->m_eqdistributions[dir] = SubDomain_.m_velSet->GetWeight(dir) *
                    SubDomain_.m_lattice[i][j][k]->density * (1.0 + 3.0 * du + 4.5 * du * du - 1.5 * u_sqr);
            }

            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[i][j][k]->m_newDistributions[dir] = SubDomain_.m_lattice[i][j][k]->m_distributions[dir] +
                    (1.0 / m_relaxation) * (SubDomain_.m_lattice[i][j][k]->m_eqdistributions[dir] - SubDomain_.m_lattice[i][j][k]->m_distributions[dir]);
            }
        }
    }
}





void Simulation::Streaming()
{
    for (int f = 0; f < SubDomain_.Involved_Lat_.size(); f++) {

        int i = SubDomain_.Involved_Lat_[f]->x_position;
        int j = SubDomain_.Involved_Lat_[f]->y_position;
        int k = SubDomain_.Involved_Lat_[f]->z_position;
        for (int dir = 0; dir < _nLatNodes; dir++) {

            if (SubDomain_.m_lattice[i][j][k]->m_neighbours[dir] != nullptr &&
                SubDomain_.m_lattice[i][j][k]->m_neighbours[dir]->NodeType_ != NodeType::Halo &&
                SubDomain_.m_lattice[i][j][k]->m_newDistributions[dir] != 0.0) {
                SubDomain_.m_lattice[i][j][k]->Stream(dir);
            }
        }
    }

    //update halo

    for (int f = 0; f < SubDomain_.Involved_Lat_.size(); f++) {

        if (SubDomain_.Involved_Lat_[f]->NodeType_ == NodeType::Halo) {
            SubDomain_.Involved_Lat_[f]->m_distributions = SubDomain_.Involved_Lat_[f]->m_newDistributions;
        }

    }
}



void Simulation::FaceCommunicate(int iter) {


    int tag_xp = 0;
    int tag_xm = 1;
    int tag_yp = 2;
    int tag_ym = 3;
    int tag_zp = 4;
    int tag_zm = 5;

    int mpi_error;
    MPI_Status status;

    

    if (SubDomain_.my_left_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                if (SubDomain_.m_lattice[1][j][k] == nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(0.0);
                    }
                }
                else {

                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(SubDomain_.m_lattice[1][j][k]->m_newDistributions[dir]);
                    }
                }
            }
        }

        

        std::vector<double> recv_data(send_data.size());
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_left_, tag_xm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_left_, tag_xp,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }

        

        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.m_lattice[0][j][k]->m_newDistributions[dir] = recv_data[index];
                    index++;
                }
            }
        }
    }

    if (SubDomain_.my_right_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                if (SubDomain_.m_lattice[SubDomain_.my_Nx_ - 2][j][k] == nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(0.0);
                    }
                }
                else {

                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(SubDomain_.m_lattice[SubDomain_.my_Nx_ - 2][j][k]->m_newDistributions[dir]);
                    }
                }
            }
        }

        

        std::vector<double> recv_data(send_data.size());
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_right_, tag_xp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_right_, tag_xm,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }

        

        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.m_lattice[SubDomain_.my_Nx_ - 1][j][k]->m_newDistributions[dir] = recv_data[index];
                    index++;
                }
            }
        }
    }


    if (SubDomain_.my_front_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                if (SubDomain_.m_lattice[i][1][k] == nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(0.0);
                    }
                }
                else {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(SubDomain_.m_lattice[i][1][k]->m_newDistributions[dir]);
                    }
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_front_, tag_ym,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_front_, tag_yp,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }

        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.m_lattice[i][0][k]->m_newDistributions[dir] = recv_data[index];
                    index++;
                }
            }
        }
    }


    if (SubDomain_.my_back_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {

                if (SubDomain_.m_lattice[i][SubDomain_.my_Ny_ - 2][k] == nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(0.0);
                    }
                }
                else {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(SubDomain_.m_lattice[i][SubDomain_.my_Ny_ - 2][k]->m_newDistributions[dir]);
                    }
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_back_, tag_yp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_back_, tag_ym,
            MPI_COMM_WORLD, &status);
        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.m_lattice[i][SubDomain_.my_Ny_ - 1][k]->m_newDistributions[dir] = recv_data[index];
                    index++;
                }
            }
        }
    }


    if (SubDomain_.my_up_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {

                if (SubDomain_.m_lattice[i][j][SubDomain_.my_Nz_ - 2] == nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(0.0);
                    }
                }
                else {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(SubDomain_.m_lattice[i][j][SubDomain_.my_Nz_ - 2]->m_newDistributions[dir]);
                    }
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_up_, tag_zp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_up_, tag_zm,
            MPI_COMM_WORLD, &status);
        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }
        int index = 0;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.m_lattice[i][j][SubDomain_.my_Nz_ - 1]->m_newDistributions[dir] = recv_data[index];
                    index++;
                }
            }
        }
    }


    if (SubDomain_.my_down_ != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {

                if (SubDomain_.m_lattice[i][j][1] == nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(0.0);
                    }
                }

                else {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(SubDomain_.m_lattice[i][j][1]->m_newDistributions[dir]);
                    }
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_down_, tag_zm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_down_, tag_zp,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }
        int index = 0;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.m_lattice[i][j][0]->m_newDistributions[dir] = recv_data[index];
                    index++;
                }
            }
        }
    }
}


double Simulation::getInletVelocity(double x, double period) {
    double vel = 0.0;
    double normalizedX = fmod(x, period);

    if (normalizedX >= 0 && normalizedX < (0.1 * period)) {
        vel = Pulsatile_.InletVelequation1(normalizedX);
    }
    else if (normalizedX >= (0.1 * period) && normalizedX < (0.3 * period)) {
        vel = Pulsatile_.InletVelequation2(normalizedX);
    }
    else if (normalizedX >= (0.3 * period) && normalizedX < (0.9 * period)) {
        vel = Pulsatile_.InletVelequation3(normalizedX);
    }

    return vel;
}





void Simulation::Iterate()
{
    total_runtime = 0.0;
    int iter = 0;
    const double tolerance = m_dt / 2; 

    while (total_runtime <= 2.7) {

        Collision();
        FaceCommunicate(iter);
        Streaming();
        InletVelocity_ = getInletVelocity(total_runtime, 0.9);
       // InletVelocity_ = 0.1;
        ApplyZouHeBc();
        ApplyOutletBc();
        ApplyWallBc();
        ComputeWallShearStress(iter);
        ComputeOSIAndTAWSS(iter);

//        if (fabs(total_runtime - 0.15) < tolerance || fabs(total_runtime - 0.9) < tolerance) {
        if (iter%1000==0 && iter!=0) {
            WriteOutput(iter);
        }

        if (myrank_ == 0) {
            std::cout << std::endl;
            std::cout << "Iteration " << iter << " has been performed successfully" << std::endl;
        }

        total_runtime += m_dt;
        iter++;
    }
}



void Simulation::ApplyZouHeBc() {

    
    for (int lat = 0; lat < SubDomain_.InletBoundaryNode_.size(); lat++) {

        SubDomain_.InletBoundaryNode_[lat]->ApplyBc(InletVelocity_);

    }
}



void Simulation::ApplyOutletBc() {

    for (int lat = 0; lat < SubDomain_.OutletBoundaryNode1_.size(); lat++) {

        for (int dir = 0; dir < _nLatNodes; dir++) {
            int x_pos = SubDomain_.OutletBoundaryNode1_[lat]->x_position;
            int y_pos = SubDomain_.OutletBoundaryNode1_[lat]->y_position;
            int z_pos = SubDomain_.OutletBoundaryNode1_[lat]->z_position;
            int opp_dir = SubDomain_.GetVelSet()->OppositeDirection(dir);
            if (SubDomain_.m_lattice[x_pos-1][y_pos][z_pos]!=nullptr){

            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[opp_dir] =
                SubDomain_.m_lattice[x_pos-1][y_pos][z_pos]->m_distributions[opp_dir];
           }
        }

    }

    for (int lat = 0; lat < SubDomain_.OutletBoundaryNode2_.size(); lat++) {

        for (int dir = 0; dir < _nLatNodes; dir++) {
            int x_pos = SubDomain_.OutletBoundaryNode2_[lat]->x_position;
            int y_pos = SubDomain_.OutletBoundaryNode2_[lat]->y_position;
            int z_pos = SubDomain_.OutletBoundaryNode2_[lat]->z_position;
            int opp_dir = SubDomain_.GetVelSet()->OppositeDirection(dir);
            if (SubDomain_.m_lattice[x_pos+1][y_pos][z_pos]!=nullptr){

            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[opp_dir] =
                SubDomain_.m_lattice[x_pos+1][y_pos][z_pos]->m_distributions[opp_dir];

          }
       }
    }

    for (int lat = 0; lat < SubDomain_.OutletBoundaryNode3_.size(); lat++) {

        for (int dir = 0; dir < _nLatNodes; dir++) {
            int x_pos = SubDomain_.OutletBoundaryNode3_[lat]->x_position;
            int y_pos = SubDomain_.OutletBoundaryNode3_[lat]->y_position;
            int z_pos = SubDomain_.OutletBoundaryNode3_[lat]->z_position;
            int opp_dir = SubDomain_.GetVelSet()->OppositeDirection(dir);
            if (SubDomain_.m_lattice[x_pos][y_pos-1][z_pos]!=nullptr){

            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[opp_dir] =
                SubDomain_.m_lattice[x_pos][y_pos-1][z_pos]->m_distributions[opp_dir];

          }
       }
    }


    for (int lat = 0; lat < SubDomain_.OutletBoundaryNode4_.size(); lat++) {

        for (int dir = 0; dir < _nLatNodes; dir++) {
            int x_pos = SubDomain_.OutletBoundaryNode4_[lat]->x_position;
            int y_pos = SubDomain_.OutletBoundaryNode4_[lat]->y_position;
            int z_pos = SubDomain_.OutletBoundaryNode4_[lat]->z_position;
            int opp_dir = SubDomain_.GetVelSet()->OppositeDirection(dir);
            if (SubDomain_.m_lattice[x_pos][y_pos-1][z_pos]!=nullptr){

            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[opp_dir] =
                SubDomain_.m_lattice[x_pos][y_pos-1][z_pos]->m_distributions[opp_dir];

          }
       }
    }
}



void Simulation::ApplyWallBc() {

    for (int lat = 0; lat < SubDomain_.WallBoundaryNode_.size(); lat++) {

        for (int dir = 0; dir < _nLatNodes; dir++) {

            if (SubDomain_.WallBoundaryNode_[lat]->m_neighbours[dir] == nullptr ||
                (SubDomain_.WallBoundaryNode_[lat]->m_neighbours[dir]->NodeType_ == NodeType::Halo &&
                    SubDomain_.WallBoundaryNode_[lat]->m_neighbours[dir]->density == 0)) {

                int opp_dir = SubDomain_.GetVelSet()->OppositeDirection(dir);

                SubDomain_.WallBoundaryNode_[lat]->m_distributions[opp_dir] =
                    SubDomain_.WallBoundaryNode_[lat]->m_distributions[dir];
            }
        }
    }
}



void Simulation::ComputeWallShearStress(int iteration) {
    double phys_den = 1060;
    double dtdt = m_dt * m_dt;
    double dxdx = SubDomain_.Geometry_.spacing * SubDomain_.Geometry_.spacing * 1e-6;
    double C = phys_den * dtdt / dxdx;
    int dim = 3;
    for (int lat = 0; lat < SubDomain_.WallBoundaryNode_.size(); lat++) {
        auto& node = SubDomain_.WallBoundaryNode_[lat];
        std::array<double, 3> tau_i = { 0.0, 0.0, 0.0 };
        for (int i = 0; i < dim; i++) {
            double sum_alpha = 0.0;
            for (int alpha = 0; alpha < _nLatNodes; alpha++) {
                double sum_k = 0.0;
                for (int k = 0; k < dim; k++) {
                    sum_k += SubDomain_.GetVelSet()->GetDirection(alpha)[k] * node->normal[k] * node->normal[i];
                }
                double sum_j = 0.0;
                for (int j = 0; j < dim; j++) {
                    sum_j += SubDomain_.GetVelSet()->GetDirection(alpha)[j] * node->normal[j] * (SubDomain_.GetVelSet()->GetDirection(alpha)[i] - sum_k);
                }
                double f_ne_alpha = node->m_distributions[alpha] - node->m_eqdistributions[alpha]; 
                       sum_alpha += f_ne_alpha * sum_j;
            }
            tau_i[i] = C * sum_alpha;
        }
        double WSS = sqrt(tau_i[0] * tau_i[0] + tau_i[1] * tau_i[1] + tau_i[2] * tau_i[2]); node->wall_shear_stress = WSS;
        for (int i = 0; i < dim; i++) {
            node->accumulated_shear_stress[i] += tau_i[i];
        }
        node->accumulated_shear_stress_magnitude += WSS;
    }
}




void Simulation::ComputeOSIAndTAWSS(int iteration) {

    for (auto& node : SubDomain_.WallBoundaryNode_) {

        double tawss = node->accumulated_shear_stress_magnitude / Pulsatile_.period;

        double norm_accumulated_shear_stress = sqrt(pow(node->accumulated_shear_stress[0], 2) +
                                                    pow(node->accumulated_shear_stress[1], 2) +
                                                    pow(node->accumulated_shear_stress[2], 2));


        double osi = 0.5 * (1.0 - (norm_accumulated_shear_stress / node->accumulated_shear_stress_magnitude));

        node->osi = osi;
        node->tawss = tawss*m_dt;
    }
}



void Simulation::SetRelaxation(double tau) {
      m_relaxation = tau;
}



void Simulation::WriteOutput(int iter) {

    OutputFile file(iter);
    file.WriteVTK_Collective(SubDomain_);
    file.WriteSurfaceAttribute(SubDomain_, iter);
    std::string inputFileName = "FlowHemodynamics/" + std::to_string(iter) + "_FlowHemodynamics.csv";
    std::string outputFolder = "FlowHemodynamics_Smoothed";
    std::string outputFileName = outputFolder + "/" + std::to_string(iter) + "_FlowHemodynamics_Smoothed.csv";
    double neighborSearchRadius = 0.5;
    file.processFlowData(inputFileName, outputFileName, neighborSearchRadius);

}

