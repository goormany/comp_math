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

// f(x) = ln(ln x) - exp(-x^2)
double f(double x) {
    return log(log(x)) - exp(-(x * x));
}


// phi(x) = x - 2.754*( ln(ln x) - exp(-x^2) )
double phi(double x) {
    if (x <= 1.0) {
        cerr << "Ошибка: x <= 1 в phi(x)" << endl;
        return x;
    }
    return x - 2.754 * (log(log(x)) - exp(-(x * x)));
}

//  phi'(x) = 1 - 2.754/(x ln x) - 5.508*x*exp(-x^2)
double phi_prime(double x) {
    if (x <= 1.0) {
        cerr << "Ошибка: x <= 1 в phi_prime(x)" << endl;
        return 0.0;
    }
    double ex2 = exp(x * x);         
    double e_neg_x2 = exp(-x * x);   
    double lnx = log(x);
    double numerator = e_neg_x2 * (2754.0 * x * x * lnx + 1377.0 * ex2);
    double denominator = 500.0 * x * lnx;
    return 1.0 - numerator / denominator;
}

double f_round(double x, double delta) {
    return Round(f(x), delta);
}

double phi_round(double x, double delta) {
    return Round(phi(x), delta);
}

double phi_prime_round(double x, double delta) {
    return Round(phi_prime(x), delta);
}


double SimpleIteration(double x0, double eps, double delta, int &iter) {
    const int MAX_ITER = 1000000;
    iter = 0;
    double x = x0;
    double x_prev;

    do {
        x_prev = x;
        double ph;
        if (delta > 0.0) {
            ph = phi_round(x, delta);
        } else {
            ph = phi(x);
        }
        x = ph;
        iter++;

        if (x <= 1.0) {
            cout << "Внимание: значение x стало <= 1 на итерации " << iter << endl;
            break;
        }

        if (fabs(x - x_prev) < eps) {
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

    double left = 2.7;
    double right = 2.8;
    double x0 = 2.7;         
    double eps = 1e-6;         
    double delta = 1e-8;      

    cout << "Метод простых итераций (релаксации) для уравнения f(x) = ln(ln x) - exp(-x^2) = 0\n";
    cout << "Отрезок, содержащий корень: [" << left << ", " << right << "]\n";
    cout << "Итерационная функция phi(x) = x - 2.754*( ln(ln x) - exp(-x^2) )\n";
    cout << "Рекомендуемое начальное приближение x0 = " << x0 << "\n\n";


    double q_max = 0.0;
    for (double x = left; x <= right; x += 0.01) {
        double q = fabs(phi_prime(x));
        if (q > q_max) q_max = q;
    }

    cout << "Максимальное |phi'(x)| на отрезке: " << q_max << " (ожидается ≈ 0.04, должно быть < 1)\n\n";

    int choice;
    cout << "Выберите режим работы:\n";
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
            double root = SimpleIteration(x0, eps, 0.0, iter);
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
            double root = SimpleIteration(x0, eps, delta, iter);
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
                double root = SimpleIteration(x0, e, 0.0, iter);
                cout << "eps = " << e << ": итераций = " << iter
                     << ", корень = " << root << endl;
            }
            break;
        }

        case 4: {
            cout << "Исследование обусловленности (влияние точности округления)\n";
            cout << "Используется eps = " << eps << endl;
            double delta_vals[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
            for (double d : delta_vals) {
                int iter;
                double root = SimpleIteration(x0, eps, d, iter);
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