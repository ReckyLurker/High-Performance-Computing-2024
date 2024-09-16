/*
    Functions to load a vector from a txt file in ASCII format 
*/

#include <vector> 
#include <iostream> 
#include <fstream>
#include "loadVector.hpp"

std::vector<double> loadVector(const char* filename){
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::string error = std::string(RED) + std::string(filename) + std::string(" couldn't be opened") + std::string(NC) + std::string("\n"); 
        throw std::runtime_error(error);
    }
 
    std::string line;
    std::vector<double> row; 
    double value; 

    while(std::getline(file, line)){
        std::istringstream iss(line);
        while(iss >> value){
            row.push_back(value);
        }
    }
    return row; // The source vector 
}
