#ifndef LBM_H
#define LBM_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>

std::map<std::string, std::string> readInputFile(const std::string& fileName) {
    std::map<std::string, std::string> inputData;
    std::ifstream inputFile(fileName);

    if (inputFile.is_open()) {
        std::string line;
        while (std::getline(inputFile, line)) {
            std::istringstream iss(line);
            std::string key;
            std::string value;

            if (std::getline(iss, key, '=')) {
                key.erase(key.find_last_not_of(" \t") + 1);

                if (std::getline(iss, value)) {
                    value.erase(0, value.find_first_not_of(" \t"));
                    value.erase(value.find_last_not_of(" \t") + 1);

                    inputData[key] = value;
                }
            }
        }
        inputFile.close();
    } else {
        std::cerr << "Unable to open file " << fileName << std::endl;
    }

    return inputData;
}

void readAndStoreVariables(std::map<std::string, std::string>& inputData,
                           double& dt, double& rho, double& nu_blood,
                           double& tau, double& dx, double& u,
                           double& reynoldsNumber, double& pressure,int num_Out,
                           std::vector<double>& inletOutletCoords) {
    if (inputData.find("time resulation") != inputData.end()) {
        dt = std::stod(inputData["time resulation"]);
    }

    if (inputData.find("density") != inputData.end()) {
        rho = std::stod(inputData["density"]);
    }

    if (inputData.find("kinematic Viscosity") != inputData.end()) {
        nu_blood = std::stod(inputData["kinematic Viscosity"]);
    }

    if (inputData.find("relaxation factor") != inputData.end()) {
        tau = std::stod(inputData["relaxation factor"]);
    }

    if (inputData.find("grid size") != inputData.end()) {
        dx = std::stod(inputData["grid size"]);
    }

    if (inputData.find("Inlet Velocity") != inputData.end()) {
        u = std::stod(inputData["Inlet Velocity"]);
    }

    if (inputData.find("Reynolds number") != inputData.end()) {
        reynoldsNumber = std::stod(inputData["Reynolds number"]);
    }

    if (inputData.find("pressure") != inputData.end()) {
        pressure = std::stod(inputData["pressure"]);
    }

    if (inputData.find("number of outlets") != inputData.end()) {
        num_Out = std::stod(inputData["number of outlets"]);
    }

    if (inputData.find("inlet/outlet coordinates") != inputData.end()) {
        std::istringstream coordStream(inputData["inlet/outlet coordinates"]);
        std::string coord;
        while (std::getline(coordStream, coord, ',')) {
            inletOutletCoords.push_back(std::stod(coord));
        }
    }
}

void printVariables(double dt, double rho, double nu_blood, double tau, double dx, double u, double reynoldsNumber,
                    double pressure,int num_Out, const std::vector<double>& inletOutletCoords) {
    std::cout << "dt: " << dt << std::endl;
    std::cout << "rho: " << rho << std::endl;
    std::cout << "nu_blood: " << nu_blood << std::endl;
    std::cout << "tau: " << tau << std::endl;
    std::cout << "dx: " << dx << std::endl;
    std::cout << "u: " << u << std::endl;
    std::cout << "Reynolds Number: " << reynoldsNumber << std::endl;
    std::cout << "pressure: " << pressure << std::endl;
    std::cout << "number of outlets: " << num_Out << std::endl;
    std::cout << "Inlet/Outlet Coordinates: ";
    for (const double& coord : inletOutletCoords) {
        std::cout << coord << " ";
    }
    std::cout << std::endl;
}

#endif

