#include <iostream>
#include <cmath>
#include <iomanip>
#include <functional>
#include <fstream>
using namespace std;


double f(double x) {
    return cosh(x * x);
}

double rectangleMethod(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    
    for (int i = 0; i < n; i++) {
        double x_mid = a + (i + 0.5) * h;
        sum += f(x_mid);
    }
    return h * sum;
}

double trapezoidalMethod(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    
    for (int i = 0; i < n; i++) {
        double x_i = a + i * h;
        double x_i1 = a + (i + 1) * h;
        sum += f(x_i) + f(x_i1);
    }
    
    return (h / 2.0) * sum;
}

double simpsonMethod(double a, double b, int n) {
    if (n % 2 != 0) {
        n++;
    }
    
    double h = (b - a) / n;
    double sum = 0.0;
    
    for (int i = 0; i < n; i += 2) {
        double x0 = a + i * h;
        double x1 = a + (i + 1) * h;
        double x2 = a + (i + 2) * h;
        sum += f(x0) + 4.0 * f(x1) + f(x2);
    }
    
    return (h / 3.0) * sum;
}

vector<int> Runge(const std::string& methodName, std::function<double(double, double, int)> method, double a, double b, double epsilon, vector<int> vals) {
    int n = 2;
    double I_h, I_h2, error;
    int iteration = 0;
    
    std::cout << "\n--- " << methodName << " ---" << std::endl;
    std::cout << std::setw(10) << "n" << std::setw(20) << "Integral" << std::setw(20) << "Error" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    
    while (true) {
        I_h = method(a, b, n);
        I_h2 = method(a, b, 2 * n);
        
        if (methodName == "Rectangle Method") {
            error = std::abs(I_h2 - I_h) / 3.0;
        } else if (methodName == "Trapezoidal Method") {
            error = std::abs(I_h2 - I_h) / 3.0;
        } else if (methodName == "Simpson Method") {
            error = std::abs(I_h2 - I_h) / 15.0;
        }
        
        std::cout << std::setw(10) << n 
                  << std::setw(20) << std::fixed << std::setprecision(10) << I_h
                  << std::setw(20) << std::scientific << std::setprecision(4) << error 
                  << std::endl;
        
        if (error < epsilon) {
            double finalIntegral;
            if (methodName == "Rectangle Method") {
                finalIntegral = I_h2 + std::abs(I_h2 - I_h) / 3.0;
            } else if (methodName == "Trapezoidal Method") {
                finalIntegral = I_h2 - std::abs(I_h2 - I_h) / 3.0;
            } else if (methodName == "Simpson Method") {
                finalIntegral = I_h2 - std::abs(I_h2 - I_h);
            }

            std::cout << std::string(50, '-') << std::endl;
            std::cout << "Result = " << std::fixed << std::setprecision(10) << finalIntegral << std::endl;
            std::cout << "n = " << n << " (h = " << (b - a) / n << ")" << std::endl;
            std::cout << "Error: " << std::scientific << std::setprecision(4) << error << std::endl;
            break;
        }
        vals.push_back(n);
        n *= 2;
    }
    return vals;
}

void to_file(string msg, vector<int> vals = {}){
    string filename = "out.txt";
    ofstream f(filename, ios::app);
    if (!f.is_open()) {
        cerr << "Ошибка: не удалось создать файл " << filename << endl;
        return;
    }
    f << msg << endl;
    for(int val: vals){
        f << val << endl;
    }
    for(int i = 0; i < 50; i++){
        f << '-';
    }
    f << endl;
    f.close();
}

int main() {
    system("rm -f out.txt");
    double a = 0.0;
    double b = 1.0;
    
    vector<double> epsilon = {1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
    for(auto eps: epsilon){
        ostringstream oss;
        oss << scientific << setprecision(0) << eps;
        string msg = "e = " + oss.str();
        to_file(msg);
        
        std::cout << "\n\n" << std::string(70, '=') << std::endl;
        std::cout << "Epsilon: ε = " << eps << std::endl;
        std::cout << std::string(70, '=') << std::endl;
        
        vector<int> vals;
        vals = Runge("Rectangle Method", rectangleMethod, a, b, eps, vals);
        to_file("Rectangle Method", vals);
        vals.clear();

        vals = Runge("Trapezoidal Method", trapezoidalMethod, a, b, eps, vals);
        to_file("Trapezoidal Method", vals);
        vals.clear();

        vals = Runge("Simpson Method", simpsonMethod, a, b, eps, vals);
        to_file("Simpson Method", vals);
        vals.clear();
    }
    
    return 0;
}