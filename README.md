# Console-Application-Development-Using-Numerical-Methods
CSE 2208 Course Assignment on Numerical Methods

## Function description

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

double rk_step(double x, double y, double h, function<double(double, double)> f): This function performs a single step of the Runge-Kutta (4th-order) method to approximate the solution of a differential equation. It takes as input the current values of x and y, the step size h and a function `f` representing the differential equation dy/dx=f(x,y). It returns the updated value of y after the step.

function<double(double, double)> equation_ax_by(double a, double b): This function returns a lambda function representing a linear differential equation of the form dy/dx=ax+by. It takes constants a and b as input and returns a lambda function that calculates the equation given values of x and y.

function<double(double, double)> equation_a_sin_b_cos(double a, double b, double c, double d): This function returns a lambda function representing a differential equation of the form dy/dx=asin(cx)+bcos(dy). It takes constants a,b,c and d as input and returns a lambda function that calculates this expression given x and y.

function<double(double, double)> equation_a_xy(double a): This function returns a lambda function for a differential equation of the form dy/dx=axy. It takes a constant a as input and returns a lambda function that evaluates this expression with given values of x and y.

function<double(double, double)> equation_a_cos_b_sin(double a, double b, double c, double d): This function returns a lambda function representing a differential equation of the form dy/dx=acos(cx)+bsin(dy). It takes constants a,b,c and d as input and returns a lambda function for the expression.

function<double(double, double)> equation_a_sin_b_y(double a, double b, double c): This function returns a lambda function for a differential equation of the form dy/dx= asin(cx)+by. It takes constants a,b and c as input and returns a lambda function for the equation with values of x and y.

function<double(double, double)> equation_a_cos_b_y(double a, double b, double c): This function returns a lambda function representing a differential equation of the form dy/dx=acos(cx)+by. It takes constants a,b and c as input and returns a lambda function to calculate the expression given x and y.

function<double(double, double)> equation_ax_b_cos(double a, double b, double d): This function returns a lambda function for a differential equation of the form dy/dx=ax+bcos(dy). It takes constants a, b and d as input and returns a lambda function for the expression.

function<double(double, double)> equation_ax_b_sin(double a, double b, double d): This function returns a lambda function representing a differential equation of the form dy/dx=ax + bsin(dy). It takes constants a,b and d as input and returns a lambda function to evaluate this expression given x and y.

function<double(double, double)> choose_equation(): This function provides a menu allowing the user to choose a differential equation type from several predefined options. It then prompts the user for any required constants based on the chosen equation type. The function returns a lambda function representing the selected differential equation, defaulting to dy/dx=ax+by if an invalid choice is entered.

void runge_kutta(): This function applies the Runge-Kutta (4th-order) method to solve a chosen differential equation. It first calls `choose_equation()` to determine the differential equation. Then it prompts the user for initial values of x,y,h(step size) and the final value of x(end of the interval). It outputs the values of x and y at each step until reaching the final x.
