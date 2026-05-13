#include <iostream>
#include <cmath>
#include <iomanip>
#include <windows.h>
#include <clocale>


using namespace std;

double f(double x) {
    return sin(x) * exp(-x * x);
}

double method_rectangle(double a, double b, int n) {
    double h = (b - a) / n;
    double S = 0.0;
    for (int i = 0; i < n; i++) {
        double x_mid = a + (i + 0.5) * h;
        S += f(x_mid);
    }
    return S * h;
}

double method_trapezoidal(double a, double b, int n) {
    double h = (b - a) / n;
    double S = (f(a) + f(b)) / 2.0;
    for (int i = 1; i < n; i++) {
        S += f(a + i * h);
    }
    return S * h;
}

double method_simpson(double a, double b, int n) {
    if (n % 2 != 0) n++;
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += (i % 2 == 0) ? 2 * f(x) : 4 * f(x);
    }
    return (h / 3.0) * sum;
}

double runge_error_rectangle(double a, double b, int n) {
    double I_n = method_rectangle(a, b, n);
    double I_2n = method_rectangle(a, b, 2 * n);
    return fabs(I_2n - I_n) / 3.0;
}

double runge_error_trapezoidal(double a, double b, int n) {
    double I_n = method_trapezoidal(a, b, n);
    double I_2n = method_trapezoidal(a, b, 2 * n);
    return fabs(I_2n - I_n) / 3.0;
}

double runge_error_simpson(double a, double b, int n) {
    double I_n = method_simpson(a, b, n);
    double I_2n = method_simpson(a, b, 2 * n);
    return fabs(I_2n - I_n) / 15.0;
}

int find_n_rectangle(double a, double b, double eps) {
    int n = 2;
    while (runge_error_rectangle(a, b, n) > eps) {
        n *= 2;
    }
    return n;
}

int find_n_trapezoidal(double a, double b, double eps) {
    int n = 2;
    while (runge_error_trapezoidal(a, b, n) > eps) {
        n *= 2;
    }
    return n;
}

int find_n_simpson(double a, double b, double eps) {
    int n = 2;
    while (runge_error_simpson(a, b, n) > eps) {
        n *= 2;
    }
    return n;
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    setlocale(LC_ALL, ".UTF-8");

    double a = 0.0, b = 1.0;
    double eps;
    
    cout << "Вычисление интеграла от sin(x) * exp(-x^2) на [0, 1]" << endl;
    cout << "Введите точность: ";
    cin >> eps;
    
    cout << endl;
    
    int n_rect = find_n_rectangle(a, b, eps);
    int n_trap = find_n_trapezoidal(a, b, eps);
    int n_simp = find_n_simpson(a, b, eps);
    
    double rect = method_rectangle(a, b, n_rect);
    double trap = method_trapezoidal(a, b, n_trap);
    double simp = method_simpson(a, b, n_simp);
    
    cout << fixed << setprecision(15);
    
    cout << "Метод прямоугольников:" << endl;
    cout << "Требуется отрезков: " << n_rect << endl;
    cout << "Значение: " << rect << endl;
    cout << "Оценка погрешности: " << runge_error_rectangle(a, b, n_rect) << endl;
    cout << endl;
    
    cout << "Метод трапеций:" << endl;
    cout << "Требуется отрезков: " << n_trap << endl;
    cout << "Значение: " << trap << endl;
    cout << "Оценка погрешности: " << runge_error_trapezoidal(a, b, n_trap) << endl;
    cout << endl;
    
    cout << "Метод Симпсона:" << endl;
    cout << "Требуется отрезков: " << n_simp << endl;
    cout << "Значение: " << simp << endl;
    cout << "Оценка погрешности: " << runge_error_simpson(a, b, n_simp) << endl;
    
    return 0;
}