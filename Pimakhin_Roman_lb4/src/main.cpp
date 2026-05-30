#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <string>

using namespace std;

double f(double x) {
    return log(x) / (x + 1.0);
}

double composite_midpoint(double a, double b, int n) {
    if (n <= 0) {
        throw invalid_argument("Number of partitions n must be positive.");
    }

    double h = (b - a) / n; // шаг разбиения (длина)
    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        double x_mid = a + (i + 0.5) * h; // середина i-го подинтервала
        sum += f(x_mid);
    }

    return h * sum;
}

double composite_trapezoid(double a, double b, int n) {
    if (n <= 0) {
        throw invalid_argument("Number of partitions n must be positive.");
    }

    double h = (b - a) / n; // h - шаг разбиения
    double sum = 0.5 * (f(a) + f(b)); // начальная сумма с половинными значениями на концах

    for (int i = 1; i < n; ++i) {
        double x = a + i * h; // x - i-я точка разбиения
        sum += f(x);
    }

    return h * sum;
}

double composite_simpson(double a, double b, int n) {
    if (n <= 0) {
        throw invalid_argument("Number of partitions n must be positive.");
    }
    if (n % 2 != 0) {
        throw invalid_argument("For Simpson's rule, n must be even.");
    }

    double h = (b - a) / n; // h - шаг разбиения
    double sum_odd = 0.0; // сумма для нечетных индексов
    double sum_even = 0.0; // сумма для четных индексов (кроме 0 и n)

    for (int i = 1; i < n; ++i) {
        double x = a + i * h; // x - i-я точка разбиения
        if (i % 2 == 0) {
            sum_even += f(x);
        } else {
            sum_odd += f(x);
        }
    }

    return (h / 3.0) * (f(a) + f(b) + 4.0 * sum_odd + 2.0 * sum_even);
}

struct IntegrationResult {
    double value;
    double error;
    int n;
    int iterations;
};

typedef double (*MethodFunc)(double, double, int); // тип функции для методов численного интегрирования

IntegrationResult adaptive_runge(double a, double b, double eps, int order, MethodFunc method) {
    if (eps <= 0.0) {
        throw invalid_argument("Accuracy eps must be positive.");
    }
    if (order <= 0) {
        throw invalid_argument("Method order must be positive.");
    }

    const double denominator = pow(2.0, order) - 1.0; // для оценки погрешности Рунге
    const double min_relative_change = 1e-2; // минимальное относительное изменение для проверки сходимости

    int n = 2;
    double I_prev = method(a, b, n); // значение интеграла для текущего числа разбиений
    int iterations = 0; // количество удвоений сетки

    while (true) {
        int n2 = 2 * n;
        double I_curr = method(a, b, n2); // значение интеграла для удвоенного числа разбиений

        double error_est = fabs(I_curr - I_prev) / denominator; // оценка погрешности Рунге
        double relative_change = fabs(I_curr - I_prev) / (fabs(I_curr) + 1e-15); // относительное изменение между текущим и предыдущим значением интеграла

        if (error_est <= eps && relative_change <= min_relative_change) {
            double corrected_value = I_curr + (I_curr - I_prev) / denominator;
            return {corrected_value, error_est, n2, iterations + 1};
        }

        n = n2;
        I_prev = I_curr;
        ++iterations;

        if (n > 1 << 26) {
            throw runtime_error("Too many partitions: possible convergence issue or requested accuracy too strict.");
        }
    }
}

void print_result(const string& title, const IntegrationResult& res, double eps) {
    cout << title << "\n";
    cout << "  integral value          = " << setprecision(15) << fixed << res.value << "\n";
    cout << "  estimated error         = " << scientific << setprecision(6) << res.error << "\n"; // оценка погрешности Рунге
    cout << "  required accuracy eps   = " << scientific << setprecision(6) << eps << "\n"; // требуемая точность
    cout << "  number of partitions n  = " << res.n << "\n"; // количество разбиений, используемых для достижения требуемой точности
    cout << "  number of grid doublings = " << res.iterations << "\n"; /// количество удвоений сетки, необходимых для достижения требуемой точности
    cout << "\n";
}

int main() {
    const double a = 1.0;
    const double b = 2.0;

    cout << "Variant 12: f(x) = ln(x) / (x + 1), integral on [1, 2]\n\n";

    while (true) {
        double eps;
        cout << "eps = ";
        if (!(cin >> eps)) {
            cout << "\nInput error. Program terminating.\n";
            return 1;
        }

        if (eps <= 0.0) {
            cout << "Terminating.\n";
            break;
        }

        try {
            IntegrationResult rect = adaptive_runge(a, b, eps, 2, composite_midpoint);
            IntegrationResult trap = adaptive_runge(a, b, eps, 2, composite_trapezoid);
            IntegrationResult simp = adaptive_runge(a, b, eps, 4, composite_simpson);

            cout << "\n==============================================\n";
            cout << "Integration interval: [1, 2]\n";
            cout << "Integrand function: f(x) = ln(x) / (x + 1)\n";
            cout << "==============================================\n\n";

            print_result("Midpoint Rule", rect, eps);
            print_result("Trapezoidal Rule", trap, eps);
            print_result("Simpson's Rule", simp, eps);

        } catch (const exception& ex) {
            cerr << "Error: " << ex.what() << "\n";
        }
    }

    return 0;
}
