#ifndef DENSEMATRIX_HPP
#define DENSEMATRIX_HPP

#include "matrix.hpp"

class DenseMatrix : public Matrix {
private:
    std::vector<std::vector<double>> data;
public:
    // Constructor //
    DenseMatrix(); 
    DenseMatrix(std::size_t rows, std::size_t cols);
    DenseMatrix(std::size_t rows, std::size_t cols, std::vector<std::vector<double>> data);
    
    // Destructor //
    ~DenseMatrix() override;
    
    // Overridden functions //
    void set(std::size_t row, std::size_t col, double value) override;
    double get(std::size_t row, std::size_t col) const override;
    void print() const override;
    void print(const char* print_string) const override;
    
    // Utlity Functions //
    void load(const char* filename); 
    std::vector<double> diag();


    // Operator Overloads // 
    double& operator()(std::size_t row, std::size_t col) override; // Setting element at (row, col)
    const double operator()(std::size_t row, std::size_t col) const override; // Getting element at (row, col)

    std::vector<double>& operator[](std::size_t row) override; // Setting row
    const std::vector<double>& operator[](std::size_t row) const override; // Accessing Row 

    // Matrix Addition // 
    DenseMatrix operator+(const DenseMatrix& B) const; 

    // Matrix Multiplication //
    DenseMatrix operator*(const DenseMatrix& B) const; // Matrix Matrix Multiplication
    std::vector<double> operator*(const std::vector<double>& B) const; // Matrix Vector Multiplication  
};

#endif // DENSEMATRIX_HPP
