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
        while(fabs(b-a)>1e-4){
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
        while(fabs(b-a)>1e-4){
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
    } while (abs(x - xo) > 0.000001);
    cout << x << ' ';
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
    } while (abs(x - xo) > 0.000001);
    cout << x1 << ' ';
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
        } while (abs(x - xo) > 0.000001);
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
int main() {

    matrix_inversion();

    return 0;
}