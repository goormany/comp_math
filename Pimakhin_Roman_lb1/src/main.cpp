#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

// функция: f(x) = (1 + cos(x)) / (3 - sin(x)) - x
double f(double x) {
    return (1.0 + cos(x)) / (3.0 - sin(x)) - x;
}

// расчет f'(x)
double f1(double x) {
    double u = 1.0 + cos(x);
    double v = 3.0 - sin(x);
    double du = -sin(x);
    double dv = -cos(x);
    
    double fraction_derivative = (du * v - u * dv) / (v * v);
    return fraction_derivative - 1.0;
}

// функция: x = phi(x)
// phi(x) = x + tau * f(x)
double phi(double x, double tau) {
    return x + tau * f(x);
}

// 1. метод бисекций
double bisectionMethod(double a, double b, double eps, int &n) {
    n = 0;
    double x;
    while ((b - a) / 2.0 > eps) {
        n++;
        x = (a + b) / 2.0;
        if (f(a) * f(x) <= 0) b = x;
        else a = x;
    }
    return (a + b) / 2.0;
}

// 2. метод хорд
double chordMethod(double a, double b, double eps, int &n) {
    n = 0;
    double x_prev = a;
    double x_curr = b;
    double x_next;
    while (true) {
        n++;
        x_next = x_curr - f(x_curr) * (x_curr - x_prev) / (f(x_curr) - f(x_prev));
        if (abs(x_next - x_curr) < eps) break;
        x_prev = x_curr;
        x_curr = x_next;
    }
    return x_next;
}

// 3. метод Ньютона
double newtonMethod(double x0, double eps, int &n) {
    n = 0;
    double x_curr = x0;
    double x_next;
    while (true) {
        n++;
        x_next = x_curr - f(x_curr) / f1(x_curr);
        if (abs(x_next - x_curr) < eps) break;
        x_curr = x_next;
    }
    return x_next;
}

// 4. метод простых итераций
double iterationMethod(double x0, double eps, int &n, double tau) {
    n = 0;
    double x_curr = x0;
    double x_next;
    while (true) {
        n++;
        x_next = phi(x_curr, tau);
        if (abs(x_next - x_curr) < eps) break;
        x_curr = x_next;
    }
    return x_next;
}

int main() {
    // интервал: [0, 1]
    double a = 0.0, b = 1.0; 
    double eps;
    int choice;

    cout << "Equation: (1 + cos(x)) / (3 - sin(x)) - x = 0 on [0, 1]" << endl;
    cout << "Enter precision (eps): ";
    cin >> eps;

    cout << "Select method:\n";
    cout << "1. Bisection Method\n";
    cout << "2. Chord Method\n";
    cout << "3. Newton's Method\n";
    cout << "4. Simple Iteration Method\n";
    cout << "Choice: ";
    cin >> choice;

    int n = 0;
    double root = 0;

    switch (choice) {
        case 1: root = bisectionMethod(a, b, eps, n); break;
        case 2: root = chordMethod(a, b, eps, n); break;
        case 3: root = newtonMethod(b, eps, n); break; 
        case 4: root = iterationMethod(b, eps, n, 0.6); break; 
        default: 
            cout << "Invalid choice!" << endl;
            return 1;
    }

    cout << fixed << setprecision(10);
    cout << "\n--- Result ---" << endl;
    cout << "Root: " << root << endl;
    cout << "Iterations: " << n << endl;
    cout << "Residual f(x): " << f(root) << endl;

    return 0;
}
