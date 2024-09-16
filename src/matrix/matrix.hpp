#ifndef MATRIX_HPP 
#define MATRIX_HPP 

/*
    batch-wise load methods for armadillo matrices. 
*/ 

#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <vector>
#include <string> 
#include <armadillo> 

#define YELLOW "\033[33m" // Yellow Terminal Color
#define RED "\033[31m" // Red Terminal Color 
#define NC "\033[0m" // Reset terminal color 



// Batch wise load for dense matrices //

// This is almost 100% slower than the supplied load function // 

// void load(const char* filename, arma::Mat<double>& matrix, const size_t batchsize){
//     std::ifstream file(filename);
//     if(!file.is_open()){
//         std::cerr << RED << "Error opening file: " << filename << NC << std::endl; 
//     }
//     size_t n_rows = 0;
//     size_t n_cols = 0;

//     std::string line;
//     std::vector<std::vector<double>> batch; 

//     while(std::getline(file, line)){
//         std::istringstream iss(line);
//         std::vector<double> row; 
//         double value; 
//         while(iss >> value){
//             row.push_back(value); // push elements in batch
//         }
//         if(n_cols == 0){
//             n_cols = row.size(); // Get number of columns 
//             matrix.insert_cols(matrix.n_cols, n_cols);
//         }
//         if(row.size() != n_cols){
//             throw std::out_of_range("\033[31mIncorrect Formatting\033[0m");
//         }
//         batch.push_back(row);
//         n_rows++;

//         // Process Batch if size(batch) > batchsize wanted
//         size_t curr_batch_size = batch.size();

//         if(curr_batch_size >= batchsize){ 
//             matrix.insert_rows(matrix.n_rows, batchsize); // Insert rows in the matrix 
//             for(size_t i = 0; i < batchsize; ++i){
//                 for(size_t j = 0; j < n_cols; ++j){
//                     matrix(matrix.n_rows - batchsize + i, j) = batch[i][j];
//                 }
//             }
//             std::vector<std::vector<double>>().swap(batch); // Clear Batch 
//         }
//     }
//     if(!batch.empty()){
//         matrix.insert_rows(matrix.n_rows, batch.size());
//         for(size_t i = 0; i < batch.size() ; ++i){
//             for(size_t j = 0; j < n_cols ; ++j){
//                 matrix(matrix.n_rows - batch.size() + i, j) = batch[i][j];
//             }
//         }
//     } // Process any remaining rows 

//     file.close(); // Close file 
// }



// Default Load is slow for sparse matrices in armadillo // 

// Use of armadillo's batch wise insertion // 
arma::SpMat<double> loadCoordAscii(const char* filename){
    std::vector<long long unsigned int> location_u; 
    std::vector<long long unsigned int> location_m;
    std::vector<double> values; 

    std::ifstream file(filename);
    long long unsigned int a,b;
    double c;

    while(file >> a >> b >> c){
        location_u.push_back(a);
        location_m.push_back(b);
        values.push_back(c);
    }

    arma::Col<long long unsigned int> lu(location_u);
    arma::Col<long long unsigned int> lm(location_m);
    
    arma::Mat<long long unsigned int> location = join_rows(lu, lm).t();

    return arma::SpMat<double>(location, arma::Col<double>(values));
}

#endif 