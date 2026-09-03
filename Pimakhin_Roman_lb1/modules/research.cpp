#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

double f(double x) {
    return (1.0 + cos(x)) / (3.0 - sin(x)) - x;
}

double f1(double x) {
    double u = 1.0 + cos(x);
    double v = 3.0 - sin(x);
    double du = -sin(x);
    double dv = -cos(x);
    return ((du * v - u * dv) / (v * v)) - 1.0;
}

double phi(double x, double tau) {
    return x + tau * f(x);
}

// Function to simulate rounding errors
double roundValue(double x, double delta) {
    if (delta <= 0.0) return x;
    return delta * round(x / delta);
}

double bisectionMethod(double a, double b, double eps, int &n, double delta) {
    n = 0;
    double x;
    while ((b - a) / 2.0 > eps) {
        n++;
        x = roundValue((a + b) / 2.0, delta);
        if (f(a) * f(x) <= 0) b = x;
        else a = x;
    }
    return x;
}

double chordMethod(double a, double b, double eps, int &n, double delta) {
    n = 0;
    double x_prev = roundValue(a, delta);
    double x_curr = roundValue(b, delta);
    double x_next;
    while (n < 1000) { // безопасный лимит
        n++;
        double f_curr = roundValue(f(x_curr), delta);
        double f_prev = roundValue(f(x_prev), delta);
        x_next = roundValue(x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev), delta);
        if (abs(x_next - x_curr) < eps) break;
        x_prev = x_curr;
        x_curr = x_next;
    }
    return x_next;
}

double newtonMethod(double x0, double eps, int &n, double delta) {
    n = 0;
    double x_curr = roundValue(x0, delta);
    double x_next;
    while (n < 1000) {
        n++;
        double f_val = roundValue(f(x_curr), delta);
        double f1_val = roundValue(f1(x_curr), delta);
        x_next = roundValue(x_curr - f_val / f1_val, delta);
        if (abs(x_next - x_curr) < eps) break;
        x_curr = x_next;
    }
    return x_next;
}

double iterationMethod(double x0, double eps, int &n, double tau, double delta) {
    n = 0;
    double x_curr = roundValue(x0, delta);
    double x_next;
    while (n = 1000 || n < 1000) { // дополнительный безопасный лимит
        n++;
        x_next = roundValue(phi(x_curr, tau), delta);
        if (abs(x_next - x_curr) < eps) break;
        x_curr = x_next;
    }
    return x_next;
}

void runResearch() {
    double a = 0.0, b = 1.0;
    double tau = 0.6;
    
    vector<double> eps_values = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12};
    vector<double> delta_values = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7};

    for (int m = 1; m <= 4; m++) {
        string filename = "method" + to_string(m) + "_research.csv";
        ofstream file(filename);
        file << "eps;iterations;root;residual\n";

        for (double eps : eps_values) {
            int n = 0;
            double r = 0;
            if (m == 1) r = bisectionMethod(a, b, eps, n, 0);
            if (m == 2) r = chordMethod(a, b, eps, n, 0);
            if (m == 3) r = newtonMethod(b, eps, n, 0);
            if (m == 4) r = iterationMethod(b, eps, n, tau, 0);
            file << scientific << setprecision(12) << eps << ";" << n << ";" << r << ";" << f(r) << "\n";
        }
        file.close();
        cout << "Generated: " << filename << endl;
    }
}

int main() {
    cout << "Starting research for Variant 12..." << endl;
    runResearch();
    cout << "Research complete. CSV files generated." << endl;
    return 0;
}
