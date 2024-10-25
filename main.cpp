#include <bits/stdc++.h>
using namespace std;

double f(vector<double> &v, double x){
    double res = 0;
    for(int i = v.size(); i > 0; --i)
        res += v[i - 1] * pow(x, v.size() - i);
    return res;
}

double df(vector<double>&v,double x){
    return (f(v, x + 0.000001) - f(v, x)) / 0.000001;
}

double interval(vector<double>&v){
    return sqrt(pow((v[1] / v[0]), 2) - 2 * (v[2] / v[0]));
}

void biSection(){
    cout << "Order of equation: " << endl;
    int n;
    cin >> n;
    vector<double> v(n + 1);
    cout << "Enter " << n + 1 << " coefficients in decreasing order of powers: " << endl;
    for(double &i:v)
        cin >> i;

    while(true){
        cout << "Enter initial guesses: " << endl;
        double a, b;
        cin >> a >> b;
        if(f(v, a) * f(v, b) > 0){
            cout << "Initial guesses are wrong. Please enter again: " << endl;
            continue;
        }
        double x = 0, xo = 1e9;
        int cnt = 0;
        while(fabs(b-a)>1e-4 && cnt < 1e6){
            ++cnt;
            xo = x;
            x = (a + b) / 2;
            if(f(v, x) == 0)
                break;
            if(f(v, x) * f(v, a) < 0)
                b = x;
            else
                a = x;
        }
        if(cnt == 1e6){
            cout << "No solution found" << endl;
            cout << "Do you want to continue? (y/n): " << endl;
            char c;
            cin >> c;
            if(c == 'n')
                break;
            else if(c!='y'){
                cout << "Invalid input. Exiting..." << endl;
                break;
            }
        }
        cout << "solution: " << x << endl;
        cout << "iteration: " << cnt << endl;
        cout << "Do you want to continue? (y/n): " << endl;
        char c;
        cin >> c;
        if(c == 'n')
            break;
        else if(c!='y'){
            cout << "Invalid input. Exiting..." << endl;
            break;
        }
    }
}

void falsePostion(){
    cout << "Order of equation: " << endl;
    int n;
    cin >> n;
    vector<double> v(n + 1);
    cout << "Enter " << n + 1 << " coefficients in decreasing order of powers: " << endl;
    for(double &i:v)
        cin >> i;
    while(true){
        cout << "Enter initial guesses: " << endl;
        double a, b;
        cin >> a >> b;
        if(f(v, a) * f(v, b) > 0){
            cout << "Initial guesses are wrong. Please enter again: " << endl;
            continue;
        }
        double x = 0, xo = 1e9;
        int cnt = 0;
        while(fabs(b-a)>1e-4 && cnt < 1e6){
            ++cnt;
            xo = x;
            x = (a * f(v, b) - b * f(v, a)) / (f(v, b) - f(v, a));
            if(f(v, x) == 0)
                break;
            if(f(v, x) * f(v, a) < 0)
                b = x;
            else
                a = x;
        }
        if(cnt == 1e6){
            cout << "No solution found" << endl;
            cout << "Do you want to continue? (y/n): " << endl;
            char c;
            cin >> c;
            if(c == 'n')
                break;
            else if(c!='y'){
                cout << "Invalid input. Exiting..." << endl;
                break;
            }
        }
        cout << "solution: " << x << endl;
        cout << "iteration: " << cnt << endl;
        cout << "Do you want to continue? (y/n): " << endl;
        char c;
        cin >> c;
        if(c == 'n')
            break;
        else if(c!='y'){
            cout << "Invalid input. Exiting..." << endl;
            break;
        }
    }
}

void newtonRaphson(){
    cout << "Order of equation: " << endl;
    int n;
    cin >> n;
    vector<double> v(n + 1);
    cout << "Enter " << n + 1 << " coefficients in decreasing order of powers: " << endl;
    for(double &i:v)
        cin >> i;
    double xmax = interval(v);
    double x1 = -xmax, x = xmax, xo;
    int cnt = 0;
    do{
        ++cnt;
        xo = x;
        if(f(v,xo)==0)
            break;
        if(df(v,xo)==0){
            cout << "Error: Division by zero" << endl;
            break;
        }
        x = xo - f(v, xo) / df(v, xo);
    } while (abs(x - xo) > 0.000001 && cnt < 1e6);
    cout << x << ' ';
    cnt = 0;
    do{
        ++cnt;
        xo = x1;
        if(f(v,xo)==0)
            break;
        if(df(v,xo)==0){
            cout << "Error: Division by zero" << endl;
            break;
        }
        x1 = xo - f(v, xo) / df(v, xo);
    } while (abs(x - xo) > 0.000001 && cnt < 1e6);
    cout << x1 << ' ';
    cnt = 0;
    if(n>=3){
        double x2 = (x + x1) / 2;
        do{
            ++cnt;
            xo = x2;
            if(f(v,xo)==0)
                break;
            if(df(v,xo)==0){
                cout << "Error: Division by zero" << endl;
                break;
            }
            x2 = xo - f(v, xo) / df(v, xo);
        } while (abs(x - xo) > 0.000001 && cnt < 0);
        cout << x2 << ' ';
    }
    cout << endl;
}

void secant(){
    cout << "Order of equation: " << endl;
    int n;
    cin >> n;
    vector<double> v(n + 1);
    cout << "Enter " << n + 1 << " coefficients in decreasing order of powers: " << endl;
    for(double &i:v)
        cin >> i;
    while(true){
        cout << "Enter initial guesses: " << endl;
        double xo, x1;
        cin >> xo >> x1;

        double x2;
        int cnt = 0;
        for (int i = 0; i < 1e6;++i){
            ++cnt;
            x2 = (xo * f(v, x1) - x1 * f(v, xo)) / (f(v, x1) - f(v, xo));
            if(abs((x2 - x1)) < 0.0001)
                break;
            xo = x1;
            x1 = x2;
        }
        if(abs((x2 - x1)) >= 0.0001){
            cout << "No solution found" << endl;
            cout << "Do you want to continue? (y/n): " << endl;
            char c;
            cin >> c;
            if(c == 'n')
                break;
            else if(c!='y'){
                cout << "Invalid input. Exiting..." << endl;
                break;
            }
        }
        else{
            cout << "solution: " << x2 << endl;
            cout << "iteration: " << cnt << endl;
        }
    }
}

//Functions related to Matrix Inversion starts from here

bool gaussJordanInversion(vector<vector<double>>& matrix,vector<vector<double>>& inverse){
    int n=matrix.size();
    inverse=vector<vector<double>>(n,vector<double>(n,0));
    int i,j;
    for(i=0;i<n;i++){
        inverse[i][i]=1.0;
    }
    for(i=0;i<n;i++){
        if(matrix[i][i]==0.0){
            bool swapped=false;
            for(j=i+1;j<n;j++){
                if(matrix[j][i]!=0.0){
                    swap(matrix[i],matrix[j]);
                    swap(inverse[i],inverse[j]);
                    swapped=true;
                    break;
                }
            }
            if(!swapped)
                return false;
        }
        double diag=matrix[i][i];
        for(j=0;j<n;j++){
            matrix[i][j]/=diag;
            inverse[i][j]/=diag;
        }
        for(j=0;j<n;j++){
            if(i!=j){
                double factor=matrix[j][i];
                for(int k=0;k<n;k++){
                    matrix[j][k]-=factor*matrix[i][k];
                    inverse[j][k]-=factor*inverse[i][k];
                }
            }
        }

    }
    return true;
}

void print_matrix(vector<vector<double>>& matrix){
    int i,j,n=matrix.size();
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            cout<<setw(10)<<matrix[i][j]<<" ";
        }
        cout<<endl;
    }
}

void matrix_inversion(){
    int n,i,j;
    cout<<"Enter the size of the matrix: ";
    cin>>n;
    vector<vector<double>> matrix(n,vector<double>(n));
    vector<vector<double>> inverse;
    cout<<"Enter the elements of the matrix row wise: "<<endl;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            cin>>matrix[i][j];
        }
    }
    if(gaussJordanInversion(matrix,inverse)){
        cout<<"The inverse of the matrix is: "<<endl;
        print_matrix(inverse);
    }
    else
        cout<<"The matrix is singular and cannot be inverted."<<endl;
}

//Functions of Matrix Inversion ends here

//Functions related to Runge-Kutta Method starts from here

double rk_step(double x,double y,double h,function<double(double,double)> f){
    double k1=h*f(x,y);
    double k2=h*f(x+0.5*h,y+0.5*k1);
    double k3=h*f(x+0.5*h,y+0.5*k2);
    double k4=h*f(x+h,y+k3);
    return (y+(1.0/6.0)*(k1+2*k2+2*k3+k4));
}

function<double(double,double)> equation_ax_by(double a,double b){
    return [a,b] (double x,double y){
        return (a*x+b*y);
    };
}

function<double(double,double)> equation_a_sin_b_cos(double a,double b,double c,double d){
    return [a,b,c,d](double x,double y){
        return (a*sin(c*x)+b*cos(d*y));
    };
}

function<double(double,double)> equation_a_xy(double a){
    return [a](double x,double y){
        return (a*x*y);
    };
}

function<double(double,double)> equation_a_cos_b_sin(double a,double b,double c,double d){
    return [a,b,c,d](double x,double y){
        return (a*cos(c*x)+b*sin(d*y));
    };
}

function<double(double,double)> equation_a_sin_b_y(double a,double b,double c){
    return [a,b,c](double x,double y){
        return (a*sin(c*x)+b*y);
    };
}

function<double(double,double)> equation_a_cos_b_y(double a,double b,double c){
    return [a,b,c](double x,double y){
        return (a*cos(c*x)+b*y);
    };
}

function<double(double,double)> equation_ax_b_cos(double a,double b,double d){
    return [a,b,d](double x,double y){
        return (a*x+b*cos(d*y));;
    };
}

function<double(double,double)> equation_ax_b_sin(double a,double b,double d){
    return [a,b,d](double x,double y){
        return (a*x+b*sin(d*y));;
    };
}

function<double(double,double)>choose_equation(){
    int choice;
    cout<<"Choose a differential equation type: "<<endl;
    cout<<"1) dy/dx=ax+by"<<endl;
    cout<<"2) dy/dx=asin(px)+bcos(qy)"<<endl;
    cout<<"3) dy/dx=axy"<<endl;
    cout<<"4) dy/dx=acos(px)+bsin(qy)"<<endl;
    cout<<"5) dy/dx=asin(px)+by"<<endl;
    cout<<"6) dy/dx=acos(px)+by"<<endl;
    cout<<"7) dy/dx=ax+bcos(py)"<<endl;
    cout<<"8) dy/dx=ax+bsin(py)"<<endl;
    cout<<"Enter your choice: ";
    cin>>choice;
    switch(choice){
        case 1:
        {
            double a,b;
            cout<<"Enter a: ";
            cin>>a;
            cout<<"Enter b: ";
            cin>>b;
            return equation_ax_by(a,b);
        }
        case 2:
        {
            double a,b,p,q;
            cout<<"Enter a: ";
            cin>>a;
            cout<<"Enter b: ";
            cin>>b;
            cout<<"Enter p: ";
            cin>>p;
            cout<<"Enter q: ";
            cin>>q;
            return equation_a_sin_b_cos(a,b,p,q);
        }
        case 3:
        {
            double a;
            cout<<"Enter a: ";
            cin>>a;
            return equation_a_xy(a);
        }
        case 4:
        {
            double a,b,p,q;
            cout<<"Enter a: ";
            cin>>a;
            cout<<"Enter b: ";
            cin>>b;
            cout<<"Enter p: ";
            cin>>p;
            cout<<"Enter q: ";
            cin>>q;
            return equation_a_cos_b_sin(a,b,p,q);
        }
        case 5:
        {
            double a,b,p;
            cout<<"Enter a: ";
            cin>>a;
            cout<<"Enter b: ";
            cin>>b;
            cout<<"Enter p: ";
            cin>>p;
            return equation_a_sin_b_y(a,b,p);
        }
        case 6:
        {
            double a,b,p;
            cout<<"Enter a: ";
            cin>>a;
            cout<<"Enter b: ";
            cin>>b;
            cout<<"Enter p: ";
            cin>>p;
            return equation_a_cos_b_y(a,b,p);
        }
        case 7:
        {
             double a,b,p;
            cout<<"Enter a: ";
            cin>>a;
            cout<<"Enter b: ";
            cin>>b;
            cout<<"Enter p: ";
            cin>>p;
            return equation_ax_b_cos(a,b,p);
        }
        case 8:
        {
             double a,b,p;
            cout<<"Enter a: ";
            cin>>a;
            cout<<"Enter b: ";
            cin>>b;
            cout<<"Enter p: ";
            cin>>p;
            return equation_ax_b_sin(a,b,p);
        }
        default:
        {
            cout<<"Invalid choice! Using dy/dx=ax+by as default."<<endl;
            return equation_ax_by(1,1);
        }
    }
}

void runge_kutta(){
    function<double(double,double)> f=choose_equation();
    double x0,y0,h,xn;
    cout<<"Enter initial value of x (x0): ";
    cin>>x0;
    cout<<"Enter initial value of y (y0): ";
    cin>>y0;
    cout<<"Enter step size (h): ";
    cin>>h;
    cout<<"Enter final value of x (xn): ";
    cin>>xn;

    double x=x0,y=y0;
    cout<<fixed<<setprecision(5);
    cout<<"x     "<<"    y"<<endl;
    cout<<x<<" "<<y<<endl;
    while(x<xn){
        double newy=rk_step(x,y,h,f);
        double newx=x+h;
        if(newx>xn){
            double extra=newx-xn;
            double adjustedh=h-extra;
            newy=rk_step(x,y,adjustedh,f);
            newx=xn;
        }
        x=newx;
        y=newy;
        if(x<xn){
            cout<<x<<" "<<y<<endl;
        }
    }
}

//Functions for Runge-Kutta Method ends here

void nonlinearEquations() {
    cout << "1. Bi-Section Method" << endl;
    cout << "2. False Position Method" << endl;
    cout << "3. Newton-Raphson Method" << endl;
    cout << "4. Secant Method\n" << endl;
    cout << "Enter your choice: ";
    int choice;
    cin >> choice;
    switch (choice) {
        case 1:
            biSection();
            break;
        case 2:
            falsePostion();
            break;
        case 3:
            newtonRaphson();
            break;
        case 4:
            secant();
            break;
        default:
            cout << "Invalid choice" << endl;
            break;
    }
}

void linearEquations() {
    
}

int main() {

    cout << "Console Application of Numerical Methods" << endl;
    cout << "1. Linear Equations" << endl;
    cout << "2. Non-Linear Equations" << endl;
    cout << "3. Differential Equations" << endl;
    cout << "4. Matrix Inversion\n" << endl;

    int choice;
    cout << "Enter your choice: ";
    cin >> choice;

    switch (choice) {
        case 1:
            linearEquations();
            break;
        case 2:
            nonlinearEquations();
            break;
        case 3:
            runge_kutta();
            break;
        case 4:
            matrix_inversion();
            break;
        default:
            cout << "Invalid choice" << endl;
            break;
    }

    return 0;
}