#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
using namespace std;

#define EPS 1e-3

double f(double x){
    return exp(x) / (1.0 + x*x);
}

double rectangleRule(double a, double b, int n){
    double h = (b-a)/n;
    double sum = 0.0;  // Fixed: use double literal

    for(int i = 0; i < n; i++){
        double x_mid = a + i * h + h/2.0;
        sum += f(x_mid);
    }
    return h * sum;
}

double trapezoidRule(double a, double b, int n){  // Fixed: n is int
    double h = (b-a)/n;  // Fixed: divide by n, not h
    double sum = f(a) + f(b);

    for(int i = 1; i < n; i++){
        double x_t = a + i*h;
        sum += 2.0 * f(x_t);
    }
    return (h/2.0) * sum;
}

double simpsonRule(double a, double b, int n){
    if (n % 2 != 0){
        n++;
    }
    double h = (b-a)/n;
    double sum = f(a) + f(b);
    
    for(int i = 1; i < n; i++){
        double x_i = a + i*h;
        if(i % 2 == 0){
            sum += 2.0 * f(x_i);
        } else {
            sum += 4.0 * f(x_i);
        }
    }
    return (h/3.0) * sum;
}

// Fixed: error passed by reference (&)
double integRange(double a, double b, double eps, int method, 
                  int &n_final, double &error){  
    int n = 2;
    double I_h, I_h2;
    int maxIters = 30;

    for(int iter = 0; iter < maxIters; iter++){
        if(method == 1){
            I_h = rectangleRule(a, b, n);
            I_h2 = rectangleRule(a, b, 2*n);
            error = abs(I_h2 - I_h) / 3.0;
        } else if(method == 2){
            // Fixed: I_h uses n, I_h2 uses 2*n
            I_h = trapezoidRule(a, b, n);
            I_h2 = trapezoidRule(a, b, 2*n);
            error = abs(I_h2 - I_h) / 3.0;
        } else {
            I_h = simpsonRule(a, b, n);
            I_h2 = simpsonRule(a, b, 2*n);
            error = abs(I_h2 - I_h) / 15.0;
        }

        if(error < eps){
            n_final = 2*n;
            return I_h2;
        }
        n *= 2;
    }
    n_final = n;
    return I_h2;
}

int main(void){
    cout << fixed << setprecision(10);
    double a = 0.0;
    double b = 1.0;

    // Fixed: negative exponents for small epsilons
    vector<double> epsilons = {1e-2, 1e-3, 1e-4};  
    
    cout << "\033[4m\033[1mFunction: \033[0m" << "e^x/(1+x^2)\n";

    for(int i = 0; i < static_cast<int>(epsilons.size()); i++){
        double eps = epsilons[i];
        cout << "\033[4m\033[1m" << "Accuracy: " << eps << "\033[0m\n";
        
        int n_rect, n_trap, n_simp;
        double err_rect, err_trap, err_simp;
        
        // Rectangle Rule
        cout << "\n\033[1mRectangle rule:\033[0m\n";
        double res_rect = integRange(a, b, eps, 1, n_rect, err_rect);
        cout << "\tFinal result: I = " << res_rect << endl;
        cout << "\tNumber of intervals: n = " << n_rect << endl;
        cout << "\tEstimated error: " << err_rect << endl;

        // Trapezoid Rule - Fixed variable names
        cout << "\n\033[1mTrapezoid rule:\033[0m\n";
        double res_trap = integRange(a, b, eps, 2, n_trap, err_trap);
        cout << "\tFinal result: I = " << res_trap << endl;
        cout << "\tNumber of intervals: n = " << n_trap << endl;  // Fixed: was n_simp
        cout << "\tEstimated error: " << err_trap << endl;         // Fixed: was err_simp

        // Simpson Rule - Fixed variable names
        cout << "\n\033[1mSimpson rule:\033[0m\n";
        double res_simp = integRange(a, b, eps, 3, n_simp, err_simp);
        cout << "\tFinal result: I = " << res_simp << endl;
        cout << "\tNumber of intervals: n = " << n_simp << endl;
        cout << "\tEstimated error: " << err_simp << endl;
        
        if (i < static_cast<int>(epsilons.size()) - 1) {
            cout << endl;
        }
    }
    return 0;
}
