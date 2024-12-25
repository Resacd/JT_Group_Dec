#include <iostream>
#include "LatticeModel.h"
#include "Simulation.h"
#include "OutputFile.h"
#include <mpi.h>
#include "Subdomain.h"
#include "Input.h"
#include "ComplexGeometry.h"
#include "disclaimer.h"
#include <chrono>
#include <thread>
#include <iomanip>
#include <omp.h>
#include "Utility.h"
#include "GeometryReader.h"
#include "Input.h"
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include "LBM.h"

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    int myrank, num_procs;
    Utility Utility;
    Reader Reader;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (myrank == 0) {
        Print_disclaimer();

        omp_set_num_threads(omp_get_max_threads());
        int numThreads = omp_get_max_threads();
        omp_set_num_threads(numThreads);
        bool stopChrono = false;
        std::thread chronoThread([&]() {
            while (!stopChrono) {
                auto now = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = now - start;

                std::cout << "\rElapsed time: " << std::fixed << std::setprecision(2)
                          << std::setw(8) << elapsed.count() << " seconds" << std::flush;

                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
        });
        chronoThread.detach();
    }

    std::map<std::string, std::string> inputData = readInputFile("input.txt");
    double dt = 0.0, rho = 0.0, nu_blood = 0.0, tau = 0.0, dx = 0.0, u = 0.0;
    double reynoldsNumber = 0.0;
    double pressure = 0.0;
    int  num_Out=0;
    std::vector<double> inletOutletCoords;

    readAndStoreVariables(inputData, dt, rho, nu_blood, tau, dx, u, reynoldsNumber, pressure,num_Out, inletOutletCoords);
   
    if (myrank == 0) {
        printVariables(dt, rho, nu_blood, tau, dx, u, reynoldsNumber, pressure,num_Out, inletOutletCoords);
    }

    std::shared_ptr<LatticeModel> velSet = std::make_shared<LatticeModel>();
    Geometry Geometry("JT2011.stl" , velSet, dx, inletOutletCoords);
    if (myrank==0){
         std::cout<<"Geometry Reconstruction is done!"<<std::endl;
    }
    SubDomain SubDomain(Geometry);
   
    if (myrank==0){
         std::cout<<"Distribution is done!"<<std::endl;
    }
    Simulation(SubDomain, tau, pressure, dt);

    MPI_Finalize();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    if (myrank == 0) {
        std::cout << "\rElapsed time: " << std::fixed << std::setprecision(2)
                  << std::setw(8) << elapsed.count() << " seconds." << std::endl;
    }

    return 0;
}

