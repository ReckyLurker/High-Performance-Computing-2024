#include <iostream> // For using input-output, cout and cin. 
#include <cstring> // For Using strcmp 
#include <fstream> // For writing computed data to files 

/* TO MAINTAIN CODE-REDABILITY */
enum Variables{
    ALPHA = 1, 
    BETA = 2, 
    GAMMA = 3, 
};

/* TAKE FLOATING POINT INPUT */
void take_input(double &var, enum Variables variable){
    std::cout << "Enter ";
    switch (variable)
    {
    case ALPHA:
        std::cout << "Alpha: ";
        break;
    case BETA:
        std::cout << "Beta: ";
        break;
    case GAMMA:
        std::cout << "Gamma: ";
        break;
    }
    std::cin >> var;
    return;
}

/* TAKE INTEGER INPUT */
void take_int_input(int &var, const char* input){
    const char *empty = "?";
    if( std::strcmp(input, empty) == 0){
        std::cout << "Enter Number of points: "; 
        std::cin >> var; 
        return; 
    }
    var = std::atoi(input);
    return; 
}

/* BOUNDARY CONDITION */
struct boundary_condition{
    double alpha; 
    double beta; 
    double gamma; 
    boundary_condition(){
        take_input(alpha, ALPHA);
        take_input(beta, BETA);
        take_input(gamma, GAMMA);
    }
};

/* GENERATE FINITE DIFFERENCE BOUNDARY */
void compute_finite_diff_matrix(const boundary_condition &boundary_condition0, const boundary_condition &boundary_condition1, 
                                        double **A, double* B, int NUM_POINTS){
    
    double dx_inv = static_cast<double>(NUM_POINTS); // compute 1/dx once to avoid division overhead; 
    double dx_inv_sqr = dx_inv*dx_inv; 

    // Boundary condition @ x = 0 // 
    A[0][0] = boundary_condition0.alpha - boundary_condition0.beta*dx_inv; 
    A[0][1] = boundary_condition0.beta*dx_inv; 
    B[0] = boundary_condition0.gamma;

    double MTWO_DX_INV_SQR = -2.0*dx_inv_sqr;

    // Rest of the points. Accessed in Row-major format because C++ is row-major // 
    for(int i = 1; i < NUM_POINTS-1; ++i){
        A[i][i-1] = dx_inv_sqr; 
        A[i][i] = MTWO_DX_INV_SQR;
        A[i][i+1] = dx_inv_sqr; 
        B[i] = 0.0;
    }
    
    // Boundary condition @ x = 1 // 
    A[NUM_POINTS-1][NUM_POINTS-1] = boundary_condition1.alpha + boundary_condition1.beta*dx_inv; 
    A[NUM_POINTS-1][NUM_POINTS-2] = -boundary_condition1.beta*dx_inv; 
    B[NUM_POINTS-1] = boundary_condition1.gamma;
    return; 
}

/* WRITE COMPUTED DATA TO FILE */
void write_to_file(double **A, double *B, int NUM_POINTS){
    std::ofstream file_for_A, file_for_B; 
    file_for_A.open("A.txt");
    file_for_B.open("B.txt");

    // Write source B // 
    for(int i = 0; i < NUM_POINTS; ++i){
        file_for_B << B[i] << "\n";
    }
    file_for_B.close(); 
    // Write Matrix A // 
    for(int i = 0; i<NUM_POINTS; ++i){
        for(int j = 0; j<NUM_POINTS; ++j){
            file_for_A << A[i][j] << " ";
        }
        file_for_A << "\n";
    }
    file_for_A.close();
    return; 
}

/* WORKFLOW LOGIC */
int main(int argc, char **argv){

    /* Some basic commands for usage guides */
    const char *help = "?";
    if( argc == 2 && std::strcmp(help, argv[1]) == 0 ){
        std::cerr << "Usage: <executable_name> <NUM_POINTS>\n";
        return 1; 
    }
    
    int NUM_POINTS = 0; 
    
    /* Take Input from terminal */
    switch(argc){
        case 1: 
            take_int_input(NUM_POINTS, help);
            break; 
        case 2:
            take_int_input(NUM_POINTS, argv[1]);
            break;
    }

    /* Take Boundary condition input */
    std::cout << "Specify Boundary condition @ x = 0: \n";
    boundary_condition boundary_condition0 = boundary_condition();
    std::cout << "Specify Boundary condition @ x = 1: \n";
    boundary_condition boundary_condition1 = boundary_condition(); 

    /* Allocate memory for finite difference matrix and source */
    /* Currently storing the full matrix A, can move on to sparse matrix later */
    double **A = new double*[NUM_POINTS];
    double *B = new double[NUM_POINTS];
    for(int i = 0; i<NUM_POINTS; ++i){
        A[i] = new double[NUM_POINTS];
    }

    /* Generate finite difference matrix with boundary conditions applied and source */
    compute_finite_diff_matrix(boundary_condition0, boundary_condition1, A, B, NUM_POINTS);
    write_to_file(A, B, NUM_POINTS);

    /* Free Up memory taken by finite difference matrix and source vector */
    for(int i = 0; i<NUM_POINTS; ++i){
        delete[] A[i]; 
    }
    delete[] A;
    delete[] B;
    return 0; 
}