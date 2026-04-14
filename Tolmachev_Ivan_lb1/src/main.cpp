#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

double f(double x) {
    return atan(x) - 1.0 / x;
}

double f1(double x) {
    return 1.0 / (1.0 + x*x) + 1.0 / (x*x);
}

double phi(double x, double tau) {
    return x - tau * f(x);
}

double roundValue(double x, double delta) {
    if (x > 0.0)
        return delta * (long((x / delta) + 0.5));
    else
        return delta * (long((x / delta) - 0.5));
}

double bisectionMethod(double left_bound, double right_bound, double eps, int &N, double delta_round = 0.0) {
    double E = fabs(eps) * 2.0;
    double f_left = f(left_bound);
    double f_right = f(right_bound);
    double x = 0.5 * (left_bound + right_bound);
    double y;

    N = 0;
    if (f_left == 0.0) return left_bound;
    if (f_right == 0.0) return right_bound;

    while ((right_bound - left_bound) >= E) {
        x = 0.5 * (right_bound + left_bound);
        y = f(x);
        if (delta_round > 0.0) y = roundValue(y, delta_round);

        if (y == 0.0) return x;
        if (y * f_left < 0.0) {
            right_bound = x;
        } else {
            left_bound = x;
            f_left = y;
        }
        N++;
    }
    return x;
}

double chordMethod(double left_bound, double right_bound, double eps, int &N, double delta_round = 0.0) {
    double f_left = f(left_bound);
    double f_right = f(right_bound);
    double x, y;

    N = 0;
    if (f_left == 0.0) return left_bound;
    if (f_right == 0.0) return right_bound;

    do {
        x = left_bound - (right_bound - left_bound) * f_left / (f_right - f_left);
        y = f(x);
        if (delta_round > 0.0) y = roundValue(y, delta_round);

        if (y == 0.0) return x;
        if (y * f_left < 0.0) {
            right_bound = x;
            f_right = y;
        } else {
            left_bound = x;
            f_left = y;
        }
        N++;
    } while (fabs(y) >= eps);

    return x;
}

double newtonMethod(double x0, double eps, int &N, double delta_round = 0.0) {
    double y, y1, dx;

    N = 0;
    do {
        y = f(x0);
        if (delta_round > 0.0) y = roundValue(y, delta_round);
        if (y == 0.0) return x0;

        y1 = f1(x0);

        dx = y / y1;
        x0 = x0 - dx;
        N++;
    } while (fabs(dx) > eps);

    return x0;
}

double iterationMethod(double x0, double eps, int &N, double tau, double delta_round = 0.0) {
    double x1 = phi(x0, tau);
    double x2 = phi(x1, tau);
    if (delta_round > 0.0) {
        x1 = roundValue(x1, delta_round);
        x2 = roundValue(x2, delta_round);
    }
    N = 2;
    while ((x1 - x2)*(x1 - x2) > fabs((2*x1 - x0 - x2) * eps)) {
        x0 = x1;
        x1 = x2;
        x2 = phi(x1, tau);
        if (delta_round > 0.0) x2 = roundValue(x2, delta_round);
        N++;
    }
    return x2;
}

int main() {
    cout << "Решение нелинейного уравнения f(x) = atan(x) - 1.0/x\n";
    cout << "Доступные методы:\n";
    cout << "  1 - Бисекция\n";
    cout << "  2 - Хорд\n";
    cout << "  3 - Ньютона\n";
    cout << "  4 - Простых итераций\n";
    cout << "  5 - Выполнить все 4 метода\n";
    cout << "\n";

    int method_choice;
    cout << "Выберите метод (1..5): ";
    cin >> method_choice;

    double left_bound = 1, right_bound = 2;
    double eps, delta_round;
    double tau;
    double x0 = 1.05;

    cout << "Введите eps (точность): ";
    cin >> eps;
    cout << "Введите delta для моделирования округления (0 - без округления): ";
    cin >> delta_round;
    
    if (method_choice == 4 || method_choice == 5) {
        cout << "Введите tau: ";
        cin >> tau;
    }

    int iterations = 0;
    double root = 0.0;

    switch (method_choice) {
        case 1:
            root = bisectionMethod(left_bound, right_bound, eps, iterations, delta_round);
            cout << "\nБисекция: x = " << root << "  итераций = " << iterations << "\n";
            break;
        case 2:
            root = chordMethod(left_bound, right_bound, eps, iterations, delta_round);
            cout << "\nХорд: x = " << root << "  итераций = " << iterations << "\n";
            break;
        case 3:
            root = newtonMethod(x0, eps, iterations, delta_round);
            cout << "\nНьютон: x = " << root << "  итераций = " << iterations << "\n";
            break;
        case 4:
            root = iterationMethod(x0, eps, iterations, tau, delta_round);
            cout << "\nИтерации: x = " << root << "  итераций = " << iterations << "\n";
            break;
        case 5: {
            cout << "\n==== Бисекция ====\n";
            root = bisectionMethod(left_bound, right_bound, eps, iterations, delta_round);
            cout << "x = " << root << "  итераций = " << iterations << "\n\n";

            cout << "==== Хорд ====\n";
            root = chordMethod(left_bound, right_bound, eps, iterations, delta_round);
            cout << "x = " << root << "  итераций = " << iterations << "\n\n";

            cout << "==== Ньютон ====\n";
            root = newtonMethod(x0, eps, iterations, delta_round);
            cout << "x = " << root << "  итераций = " << iterations << "\n\n";

            cout << "==== Итерации (tau=" << tau << ") ====\n";
            root = iterationMethod(x0, eps, iterations, tau, delta_round);
            cout << "x = " << root << "  итераций = " << iterations << "\n\n";
            break;
        }
        default:
            cout << "Неверный выбор метода\n";
            return 1;
    }

    return 0;
}