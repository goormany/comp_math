#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

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
    double left_bound = 1;
    double right_bound = 2;
    double x0 = 1.05;
    double tau = 0.8;
    vector<double> eps_values = {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6};
    vector<double> delta_values = {1e-1,1e-2,1e-3,1e-4,1e-5};

    for (int method = 1; method <= 4; ++method) {
        string fname_eps = "method" + to_string(method) + "_iterations_vs_eps.csv";
        string fname_delta = "method" + to_string(method) + "_sensitivity_delta.csv";

        ofstream ofs_eps(fname_eps.c_str());
        ofstream ofs_delta(fname_delta.c_str());
        if (!ofs_eps.is_open() || !ofs_delta.is_open()) {
            cerr << "Не удалось открыть файл для записи.\n";
            return 1;
        }

        ofs_eps << "eps;iterations;root;f_root\n";
        ofs_delta << "delta;iterations;root;f_root\n";

        for (double eps : eps_values) {
            int iterations = 0;
            double root = 0.0;
            double froot = 0.0;

            switch (method) {
                case 1:
                    root = bisectionMethod(left_bound, right_bound, eps, iterations, 0.0);
                    break;
                case 2:
                    root = chordMethod(left_bound, right_bound, eps, iterations, 0.0);
                    break;
                case 3:
                    root = newtonMethod(x0, eps, iterations, 0.0);
                    break;
                case 4:
                    root = iterationMethod(x0, eps, iterations, tau, 0.0);
                    break;
            }
            froot = f(root);
            ofs_eps << std::scientific << setprecision(6) << eps << ";" << iterations << ";" 
                    << std::setprecision(10) << root << ";" << froot << "\n";
        }

        double eps_fixed = 1e-6;
        for (double delta : delta_values) {
            int iterations = 0;
            double root = 0.0;
            double froot = 0.0;

            switch (method) {
                case 1:
                    root = bisectionMethod(left_bound, right_bound, eps_fixed, iterations, delta);
                    break;
                case 2:
                    root = chordMethod(left_bound, right_bound, eps_fixed, iterations, delta);
                    break;
                case 3:
                    root = newtonMethod(x0, eps_fixed, iterations, delta);
                    break;
                case 4:
                    root = iterationMethod(x0, eps_fixed, iterations, tau, delta);
                    break;
            }
            froot = f(root);
            ofs_delta << std::scientific << setprecision(6) << delta << ";" << iterations << ";" 
                      << std::setprecision(10) << root << ";" << froot << "\n";
        }

        ofs_eps.close();
        ofs_delta.close();
    }
    return 0;
}