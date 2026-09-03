#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <string>

using namespace std;

double f(double x) {
    return sin(x * x * x);
}

double composite_midpoint(double a, double b, int n) {
    if (n <= 0) {
        throw invalid_argument("Число разбиений n должно быть положительным.");
    }

    double h = (b - a) / n;
    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        double x_mid = a + (i + 0.5) * h;
        sum += f(x_mid);
    }

    return h * sum;
}

double composite_trapezoid(double a, double b, int n) {
    if (n <= 0) {
        throw invalid_argument("Число разбиений n должно быть положительным.");
    }

    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += f(x);
    }

    return h * sum;
}

double composite_simpson(double a, double b, int n) {
    if (n <= 0) {
        throw invalid_argument("Число разбиений n должно быть положительным.");
    }
    if (n % 2 != 0) {
        throw invalid_argument("Для формулы Симпсона n должно быть чётным.");
    }

    double h = (b - a) / n;
    double sum_odd = 0.0;
    double sum_even = 0.0;

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
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

typedef double (*MethodFunc)(double, double, int);

IntegrationResult adaptive_runge(double a, double b, double eps, int order, MethodFunc method) {
    if (eps <= 0.0) {
        throw invalid_argument("Точность eps должна быть положительной.");
    }
    if (order <= 0) {
        throw invalid_argument("Порядок метода должен быть положительным.");
    }

    const double denominator = pow(2.0, order) - 1.0;
    const double min_relative_change = 1e-2;

    int n = 2;
    double I_prev = method(a, b, n);
    int iterations = 0;

    while (true) {
        int n2 = 2 * n;
        double I_curr = method(a, b, n2);

        double error_est = fabs(I_curr - I_prev) / denominator;
        double relative_change = fabs(I_curr - I_prev) / (fabs(I_curr) + 1e-15);

        if (error_est <= eps && relative_change <= min_relative_change) {
            double corrected_value = I_curr + (I_curr - I_prev) / denominator;
            return {corrected_value, error_est, n2, iterations + 1};
        }

        n = n2;
        I_prev = I_curr;
        ++iterations;

        if (n > 1 << 26) {
            throw runtime_error("Слишком большое число разбиений: возможна проблема со сходимостью или с заданной точностью.");
        }
    }
}

void print_result(const string& title, const IntegrationResult& res, double eps) {
    cout << title << "\n";
    cout << "  значение интеграла      = " << setprecision(15) << fixed << res.value << "\n";
    cout << "  оценка погрешности      = " << scientific << setprecision(6) << res.error << "\n";
    cout << "  требуемая точность eps  = " << scientific << setprecision(6) << eps << "\n";
    cout << "  число разбиений n       = " << res.n << "\n";
    cout << "  число удвоений сетки    = " << res.iterations << "\n";
    cout << "\n";
}

int main() {
    const double a = 1.0;
    const double b = 2.0;

    cout << "Вариант 21: f(x) = sin(x^3), интеграл на [1, 2]\n";

    while (true) {
        double eps;
        cout << "eps = ";
        if (!(cin >> eps)) {
            cout << "\nОшибка ввода. Завершение программы.\n";
            return 1;
        }

        if (eps <= 0.0) {
            cout << "Завершение работы.\n";
            break;
        }

        try {
            IntegrationResult rect = adaptive_runge(a, b, eps, 2, composite_midpoint);
            IntegrationResult trap = adaptive_runge(a, b, eps, 2, composite_trapezoid);
            IntegrationResult simp = adaptive_runge(a, b, eps, 4, composite_simpson);

            cout << "\n==============================================\n";
            cout << "Интервал интегрирования: [1, 2]\n";
            cout << "Подынтегральная функция: f(x) = sin(x^3)\n";
            cout << "==============================================\n\n";

            print_result("Формула прямоугольников (средние)", rect, eps);
            print_result("Формула трапеций", trap, eps);
            print_result("Формула Симпсона", simp, eps);

        } catch (const exception& ex) {
            cerr << "Ошибка: " << ex.what() << "\n";
        }
    }

    return 0;
}