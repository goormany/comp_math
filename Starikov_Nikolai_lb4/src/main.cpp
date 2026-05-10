#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x) {
    return pow(x, 4) * exp(-x * x);
}

double rectangle(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double x_mid = a + (i + 0.5) * h;
        sum += f(x_mid);
    }
    return h * sum;
}

double trapezoid(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; ++i) {
        sum += f(a + i * h);
    }
    return h * sum;
}

double simpson(double a, double b, int n) {
    if (n % 2 != 0) n++;
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += (i % 2 == 0) ? 2.0 * f(x) : 4.0 * f(x);
    }
    return h * sum / 3.0;
}

void runge_method(double a, double b, double eps, const string& name,
                  double (*method)(double, double, int), int order) {
    int n = 2;
    int factor = (order == 2) ? 3 : 15;
    double I_h, I_h2, error;
    int max_iter = 20;

    for (int iter = 0; iter < max_iter; ++iter) {
        I_h = method(a, b, n);
        I_h2 = method(a, b, 2 * n);
        error = fabs(I_h2 - I_h) / factor;
        if (error <= eps) break;
        n *= 2;
    }

    double I_improved = I_h2 + (I_h2 - I_h) / factor;
    int n_final = 2 * n;               
    double h_final = (b - a) / n_final;

    cout << fixed << setprecision(12);
    cout << name << ":\n";
    cout << "  Количество разбиений n = " << n_final << "\n";
    cout << "  Шаг интегрирования h = " << h_final << "\n";
    cout << "  Уточнённое значение интеграла = " << I_improved << "\n";
    cout << "  Оценка погрешности = " << scientific << error << "\n\n";
}

int main() {
    double a = 0.0, b = M_PI;
    double eps;

    cout << "Введите требуемую точность epsilon (например, 0.01, 0.001, 0.0001): ";
    cin >> eps;

    cout << "\n========== Интеграл: int(0..pi) x^4 * exp(-x^2) dx ==========\n";
    cout << "Точность epsilon = " << eps << "\n\n";

    runge_method(a, b, eps, "Метод прямоугольников", rectangle, 2);
    runge_method(a, b, eps, "Метод трапеций", trapezoid, 2);
    runge_method(a, b, eps, "Метод Симпсона", simpson, 4);

    return 0;
}