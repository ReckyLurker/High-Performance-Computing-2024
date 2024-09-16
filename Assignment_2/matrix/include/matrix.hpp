#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>

class Matrix {
protected:
    std::size_t n_rows; // Number of rows 
    std::size_t n_cols; // Number of columns 
public:
    // Constructor
    Matrix(std::size_t rows, std::size_t cols) : n_rows(rows), n_cols(cols) {}
    
    // Virtual destructor
    virtual ~Matrix() = default;

    // Getters
    const std::size_t getRows() const{
        return n_rows; 
    }
    const std::size_t getCols() const{
        return n_cols; 
    }
    
    // Matrix Operations // 
    virtual void set(std::size_t row, std::size_t col, double value) = 0;
    virtual double get(std::size_t row, std::size_t col) const = 0;

    // Operator Overloads // 
    virtual double& operator()(std::size_t row, std::size_t col) = 0; // Setting element at (row, col)
    virtual const double operator()(std::size_t row, std::size_t col) const = 0; // Getting element at (row, col)

    virtual std::vector<double>& operator[](std::size_t row) = 0; // Setting row
    virtual const std::vector<double>& operator[](std::size_t row) const = 0; // Accessing Row 

    // Print to terminal // 
    virtual void print() const = 0;
    virtual void print(const char* print_string) const = 0; // Print with a print_string
};

#endif // MATRIX_HPP
