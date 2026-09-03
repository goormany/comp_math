#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

double f(double x) {
    double lnx = log(x);
    return lnx * lnx - 1.0 / x;
}

double f1(double x) {
    double lnx = log(x);
    return 2.0 * lnx / x + 1.0 / (x * x);
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
    cout << "Решение нелинейного уравнения f(x) = ln^2(x) - 1/x (вариант 19)\n";
    cout << "Доступные методы:\n";
    cout << "  1 - Бисекция\n";
    cout << "  2 - Хорд\n";
    cout << "  3 - Ньютона\n";
    cout << "  4 - Простых итераций\n";
    cout << "  5 - Выполнить все 4 метода (с сохранением CSV)\n";
    cout << "\n";

    int method_choice;
    cout << "Выберите метод (1..5): ";
    cin >> method_choice;

    double left = 1.9, right = 2.2;
    double x0 = 2.0;
    double tau = 0.0;
    if (method_choice == 4 || method_choice == 5) {
        cout << "Введите параметр tau для метода простых итераций: ";
        cin >> tau;
    }
    double eps_fixed = 1e-8;

    auto investigateEps = [&](const string& filename, int method) {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Не удалось создать файл " << filename << endl;
            return;
        }
        file << "eps;iter;root\n";

        cout << "\nИсследование зависимости от eps для метода " << method << "\n";
        cout << "   eps        iterations       root\n";
        for (double eps = 0.1; eps >= 1e-6; eps /= 10.0) {
            int iter = 0;
            double root;
            if (method == 1)
                root = bisectionMethod(left, right, eps, iter, 0.0);
            else if (method == 2)
                root = chordMethod(left, right, eps, iter, 0.0);
            else if (method == 3)
                root = newtonMethod(x0, eps, iter, 0.0);
            else if (method == 4)
                root = iterationMethod(x0, eps, iter, tau, 0.0);
            else           
                return;

            cout << " " << scientific << eps
                 << "   " << setw(4) << iter
                 << "         " << fixed << root << "\n";
            file << scientific << eps << ";" << iter << ";" << fixed << root << "\n";
        }
        file.close();
        cout << "Результаты сохранены в " << filename << endl;
    };

    auto investigateDelta = [&](const string& filename, int method) {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Не удалось создать файл " << filename << endl;
            return;
        }
        file << "delta;iter;root\n";

        double deltas[] = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001};
        int num_deltas = sizeof(deltas)/sizeof(deltas[0]);

        cout << "\nИсследование влияния delta (eps = " << eps_fixed << ") для метода " << method << "\n";
        cout << "   delta       iterations       root\n";
        for (int i = 0; i < num_deltas; ++i) {
            double delta = deltas[i];
            int iter = 0;
            double root;
            if (method == 1)
                root = bisectionMethod(left, right, eps_fixed, iter, delta);
            else if (method == 2)
                root = chordMethod(left, right, eps_fixed, iter, delta);
            else if (method == 3)
                root = newtonMethod(x0, eps_fixed, iter, delta);
            else if (method == 4)
                root = iterationMethod(x0, eps_fixed, iter, tau, delta);
            else
                return;

            cout << " " << scientific << delta
                 << "   " << setw(4) << iter
                 << "         " << fixed << root << "\n";
            file << scientific << delta << ";" << iter << ";" << fixed << root << "\n";
        }
        file.close();
        cout << "Результаты сохранены в " << filename << endl;
    };

    switch (method_choice) {
        case 1:
            investigateEps("bisection_eps.csv", 1);
            investigateDelta("bisection_delta.csv", 1);
            break;
        case 2:
            investigateEps("chord_eps.csv", 2);
            investigateDelta("chord_delta.csv", 2);
            break;
        case 3:
            investigateEps("newton_eps.csv", 3);
            investigateDelta("newton_delta.csv", 3);
            break;
        case 4:
            investigateEps("iteration_eps.csv", 4);
            investigateDelta("iteration_delta.csv", 4);
            break;
        case 5:
            investigateEps("bisection_eps.csv", 1);
            investigateDelta("bisection_delta.csv", 1);
            investigateEps("chord_eps.csv", 2);
            investigateDelta("chord_delta.csv", 2);
            investigateEps("newton_eps.csv", 3);
            investigateDelta("newton_delta.csv", 3);
            investigateEps("iteration_eps.csv", 4);
            investigateDelta("iteration_delta.csv", 4);
            break;
        default:
            cout << "Неверный выбор метода\n";
            return 1;
    }

    cout << "\nВсе исследования завершены. CSV-файлы созданы.\n";
    return 0;
}
