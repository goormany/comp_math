#include <iostream>
#include <cmath>
#include <iomanip>
#include <windows.h>
#include <clocale>

using namespace std;


double Round(double x, double delta) {
    if (delta <= 0.0) return x;
    return round(x / delta) * delta;
}

// f(x) = ln(ln(x)) - exp(-x^2)
double f(double x) {
    return log(log(x)) - exp(-(x * x));
}

// f'(x) = 1/(x*ln(x)) + 2*x*exp(-x^2)
double f_prime(double x) {
    return 1.0 / (x * log(x)) + 2.0 * x * exp(-(x * x));
}

double f_round(double x, double delta) {
    return Round(f(x), delta);
}

double f_prime_round(double x, double delta) {
    return Round(f_prime(x), delta);
}

double Newton(double x0, double eps, double delta, int &iter) {
    const int MAX_ITER = 1000000;
    iter = 0;
    double x = x0;
    double x_prev;

    do {
        x_prev = x;
        double fx, fpx;

        if (delta > 0.0) {
            fx = f_round(x, delta);
            fpx = f_prime_round(x, delta);
        } else {
            fx = f(x);
            fpx = f_prime(x);
        }

        if (fabs(fpx) < 1e-12) {
            cout << "Производная близка к нулю. Метод Ньютона расходится." << endl;
            return x;
        }

        x = x - fx / fpx;
        iter++;

        if (fabs(x - x_prev) < eps || fabs(fx) < eps) {
            break;
        }

        if (iter >= MAX_ITER) {
            cout << "Достигнуто максимальное число итераций." << endl;
            break;
        }
    } while (true);

    return x;
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    setlocale(LC_ALL, ".UTF-8");

    double a = 2.0;          
    double b = 3.0;          
    double x0 = 2.0;         
    double eps = 1e-6;        
    double delta = 1e-8;      

    int choice;
    cout << "Метод Ньютона для уравнения f(x)=ln(ln(x))-exp(-x^2)=0\n";
    cout << "Отрезок, содержащий корень: [" << a << ", " << b << "]\n";
    cout << "Рекомендуемое начальное приближение x0 = " << x0 << "\n\n";
    cout << "1. Вычисление корня без округления\n";
    cout << "2. Вычисление корня с округлением (фиксированное delta)\n";
    cout << "3. Исследование скорости сходимости (без округления, варьирование eps)\n";
    cout << "4. Исследование обусловленности (влияние точности округления delta)\n";
    cout << "0. Выход\n";
    cout << "Ваш выбор: ";
    cin >> choice;

    switch (choice) {
        case 1: {
            cout << "Введите точность eps (например, 1e-6): ";
            cin >> eps;
            int iter;
            double root = Newton(x0, eps, 0.0, iter);
            cout << fixed << setprecision(10);
            cout << "\nКорень: x = " << root << endl;
            cout << "Значение f(x) = " << f(root) << endl;
            cout << "Число итераций: " << iter << endl;
            break;
        }

        case 2: {
            cout << "Введите точность eps (например, 1e-6): ";
            cin >> eps;
            cout << "Используется delta = " << delta << endl;
            int iter;
            double root = Newton(x0, eps, delta, iter);
            cout << fixed << setprecision(10);
            cout << "\nКорень: x = " << root << endl;
            cout << "f_raw(x) = " << f(root) << endl;
            cout << "f_round(x) = " << f_round(root, delta) << endl;
            cout << "Число итераций: " << iter << endl;
            break;
        }

        case 3: {
            cout << "Исследование скорости сходимости (без округления)\n";
            double eps_vals[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
            for (double e : eps_vals) {
                int iter;
                double root = Newton(x0, e, 0.0, iter);
                cout << "eps = " << e << ": итераций = " << iter
                     << ", корень = " << root << endl;
            }
            break;
        }

        case 4: {
            cout << "Исследование обусловленности (влияние точности округления)\n";
            double delta_vals[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
            for (double d : delta_vals) {
                int iter;
                double root = Newton(x0, eps, d, iter);
                cout << "delta = " << d << ": итераций = " << iter
                     << ", корень = " << root
                     << ", f_raw(x) = " << f(root)
                     << ", f_round(x) = " << f_round(root, d) << endl;
            }
            break;
        }

        case 0:
            cout << "Выход.\n";
            break;

        default:
            cout << "Неверный выбор.\n";
    }

    return 0;
}