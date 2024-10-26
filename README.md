# Console-Application-Development-Using-Numerical-Methods
CSE 2208 Course Assignment on Numerical Methods

## Function description


### Linear Equations

void Gauss_Seidel_method(): This function calculates the solution of a system of linear equations using the Gauss-Seidel method. It first prompts the user for the number of variables in the system, then reads the coefficient matrix a and constant vector b from user input. It initializes a vector x with zeros to store the initial guess for the solution. The function then prompts the user for the maximum number of iterations and the tolerance value. It iterates through the Gauss-Seidel method, updating the solution vector x until the maximum number of iterations is reached or the solution converges within the specified tolerance. If the solution converges, the function prints the solution vector x with six decimal places. If the solution does not converge within the maximum number of iterations, a message is displayed indicating that the method failed to converge.

void Jacobi_Iterative_Method(): This function calculates the solution of a system of linear equations using the Jacobi iterative method. It first prompts the user for the number of variables in the system, then reads the coefficient matrix a and constant vector b from user input. It initializes a vector x with zeros to store the initial guess for the solution. The function then prompts the user for the maximum number of iterations and the tolerance value. It iterates through the Jacobi method, updating the solution vector x until the maximum number of iterations is reached or the solution converges within the specified tolerance. If the solution converges, the function prints the solution vector x with six decimal places. If the solution does not converge within the maximum number of iterations, a message is displayed indicating that the method failed to converge.

void LU_Factorization(): This function calculates the solution of a system of linear equations using LU factorization. It first prompts the user for the number of variables in the system, then reads the coefficient matrix a and constant vector b from user input. It decomposes the coefficient matrix a into lower and upper triangular matrices using LU factorization. It then solves the system of equations using forward and backward substitution. If the solution is successful, the function prints the solution vector with six decimal places. If the matrix is singular or the factorization fails, an error message is displayed.

solveUsingGaussianElimination(vector<vector<double>>& a, vector<double>& b): This function applies Gaussian elimination with partial pivoting to solve a system of linear equations where a is a reference to a 2D vector representing the coefficient matrix and b is a reference to a vector containing the constants on the right side of each equation. Initially, the function augments matrix a by appending vector b as an extra column to prepare it for elimination. For each column, it identifies the row with the largest pivot element and swaps rows if necessary (partial pivoting) to maximize numerical accuracy. The function then performs forward elimination to create zeros below each pivot element by subtracting scaled versions of the pivot row from the rows below. Once in upper triangular form, it proceeds with back substitution to calculate each variable's value, starting from the last row upward. If any pivot value is near zero, indicating that the matrix is singular, the function outputs a warning and returns an empty vector. Otherwise it returns a vector with the solution values for each variable in the system.

gaussianElimination(): This function serves as a high-level interface for solving a system of linear equations using Gaussian elimination. It does not take any parameters. Instead, it manages input, computation and output within its scope. This function first prompts the user for the number of variables in the system, then gathers user input for the coefficient matrix a and the constant vector b. After collecting these values, it calls solveUsingGaussianElimination with a and b as arguments to perform the calculations. If a solution vector is returned, the function prints each variable's solution, formatted to six decimal places. If the matrix is singular, a message is displayed and no solution is output. Thus, gaussianElimination orchestrates the entire process of solving linear equations from receiving input to delivering results.

solveUsingGaussJordan(vector<vector<double>>& a, vector<double>& b, int& check): This function performs the Gauss-Jordan elimination process on a system of linear equations to solve for the variables, where a is a reference to a 2D vector representing the coefficient matrix, b is a reference to a vector containing constants and check is a reference to an integer that is updated to indicate if the matrix is singular. The function begins by augmenting matrix a with vector b to create an extended matrix for the elimination process. For each column, it finds the row with the largest pivot element, then swaps rows if necessary to ensure numerical stability (partial pivoting). If any pivot element is near zero, indicating a nearly singular matrix, the function outputs an error message, sets check to 0 and exits. If the pivot is non-zero, it normalizes the pivot row by dividing each element by the pivot value. Then it performs row operations to create zeros above and below each pivot, effectively reducing the matrix to reduced row echelon form (RREF). Finally the solution is stored in vector b with each element representing a solution to one variable in the system.

gaussJordan(): This function is the primary interface for solving a system of linear equations using the Gauss-Jordan elimination method. It takes no parameters, handling input, processing and output directly. The function first prompts the user for the number of variables, then reads the coefficient matrix a and constant vector b from user input. It initializes a variable check to 1 which will indicate whether the matrix is singular. The function then calls solveUsingGaussJordan, passing a, b and check as arguments to perform the elimination. If the solution is successful (check remains 1), it outputs each solution with six decimal places, formatted as x1, x2 etc. If the matrix is singular, a message is displayed and no solution is printed. Thus gaussJordan encapsulates the complete Gauss-Jordan elimination process from input to result presentation.


### Non-Linear Equations

double f(vector<double>&v,double x): This function takes a vector which contains coefficients of a polynomial and a value of x as input. It returns the value of the polynomial at x.

double df(vector<double>&v,double x): This function takes a vector which contains coefficients of a polynomial and a value of x as input. It returns the value of the derivative of the polynomial at x.

double interval(vector<double>&v): This function takes a vector which contains coefficients of a polynomial as input. It returns the interval in which the root of the polynomial lies.

void bisection(): It calculates the root of a polynomial using the bisection method. It first takes the coefficients of the polynomial as input and then continiously asks for two initial guesses and finds the root of the polynomial. The process can be terminated by given choice.

void falsePostion(): It calculates the root of a polynomial using the false position method. It first takes the coefficients of the polynomial as input and then continiously asks for two initial guesses and finds the root of the polynomial. The process can be terminated by given choice.

void newtonRaphson(): It calculates the root of a polynomial using the Newton-Raphson method. It first takes the coefficients of the polynomial as input. Then it uses the interval function to find the interval in which the root lies. Then it asks for an initial guess and finds the root of the polynomial. The process can be terminated by given choice.

void secant(): It calculates the root of a polynomial using the secant method. It first takes the coefficients of the polynomial as input. Then it asks for two initial guesses and finds the root of the polynomial. The process can be terminated by given choice.


### Matrix Inversion

bool gaussJordanInversion(vector<vector<double>>& matrix, vector<vector<double>>& inverse): This function calculates the inverse of a matrix using the Gauss-Jordan elimination method. It takes two matrices as inputs, matrix (the matrix to be inverted) and inverse (where the result will be stored). The function first initializes inverse as an identity matrix. Then it performs row operations to transform matrix into the identity matrix while performing the same operations on inverse to produce the inverse matrix. The function returns true if the inversion is successful. Otherwise, it returns false if the matrix is singular.

void print_matrix(vector<vector<double>>& matrix): This function prints a given matrix to the console. It takes a matrix as input and displays each element formatted to align neatly in columns for readability. Each element is displayed with a specific width for consistent spacing.

void matrix_inversion(): This function facilitates user input for matrix inversion. It prompts the user to input the size of the matrix and its elements row by row. It then attempts to compute the inverse of the matrix using gaussJordanInversion(). If successful, it displays the inverse. Otherwise, it informs the user that the matrix is singular and cannot be inverted.


### Differential Equations

double rk_step(double x, double y, double h, function<double(double, double)> f): This function performs a single step of the Runge-Kutta (4th-order) method to approximate the solution of a differential equation. It takes as input the current values of x and y, the step size h and a function f representing the differential equation dy/dx=f(x,y). It returns the updated value of y after the step.

function<double(double, double)> equation_ax_by(double a, double b): This function returns a lambda function representing a linear differential equation of the form dy/dx=ax+by. It takes constants a and b as input and returns a lambda function that calculates the equation given values of x and y.

function<double(double, double)> equation_a_sin_b_cos(double a, double b, double c, double d): This function returns a lambda function representing a differential equation of the form dy/dx=asin(cx)+bcos(dy). It takes constants a,b,c and d as input and returns a lambda function that calculates this expression given x and y.

function<double(double, double)> equation_a_xy(double a): This function returns a lambda function for a differential equation of the form dy/dx=axy. It takes a constant a as input and returns a lambda function that evaluates this expression with given values of x and y.

function<double(double, double)> equation_a_cos_b_sin(double a, double b, double c, double d): This function returns a lambda function representing a differential equation of the form dy/dx=acos(cx)+bsin(dy). It takes constants a,b,c and d as input and returns a lambda function for the expression.

function<double(double, double)> equation_a_sin_b_y(double a, double b, double c): This function returns a lambda function for a differential equation of the form dy/dx= asin(cx)+by. It takes constants a,b and c as input and returns a lambda function for the equation with values of x and y.

function<double(double, double)> equation_a_cos_b_y(double a, double b, double c): This function returns a lambda function representing a differential equation of the form dy/dx=acos(cx)+by. It takes constants a,b and c as input and returns a lambda function to calculate the expression given x and y.

function<double(double, double)> equation_ax_b_cos(double a, double b, double d): This function returns a lambda function for a differential equation of the form dy/dx=ax+bcos(dy). It takes constants a, b and d as input and returns a lambda function for the expression.

function<double(double, double)> equation_ax_b_sin(double a, double b, double d): This function returns a lambda function representing a differential equation of the form dy/dx=ax + bsin(dy). It takes constants a,b and d as input and returns a lambda function to evaluate this expression given x and y.

function<double(double, double)> choose_equation(): This function provides a menu allowing the user to choose a differential equation type from several predefined options. It then prompts the user for any required constants based on the chosen equation type. The function returns a lambda function representing the selected differential equation, defaulting to dy/dx=ax+by if an invalid choice is entered.

void runge_kutta(): This function applies the Runge-Kutta (4th-order) method to solve a chosen differential equation. It first calls choose_equation() to determine the differential equation. Then it prompts the user for initial values of x,y,h(step size) and the final value of x(end of the interval). It outputs the values of x and y at each step until reaching the final x.
