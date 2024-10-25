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