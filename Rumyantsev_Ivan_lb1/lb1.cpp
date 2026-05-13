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
    return asin(2.0 * x / (1.0 + x * x)) - exp(-x * x);
}


double roundValue(double x, double delta) {
    if (delta <= 0.0) return x;
    if (x > 0.0)
        return delta * (long((x / delta) + 0.5));
    else
        return delta * (long((x / delta) - 0.5));
}


double bisectionMethod(double a, double b, double eps, int &N, double delta_round = 0.0) {
    double E = fabs(eps) * 2.0;
    double fa = f(a);
    double fb = f(b);
    double x, fx;

    N = 0;
    if (fa == 0.0) return a;
    if (fb == 0.0) return b;

    while ((b - a) >= E) {
        x = 0.5 * (a + b);
        fx = f(x);
        if (delta_round > 0.0) fx = roundValue(fx, delta_round);

        if (fx == 0.0) return x;
        if (fa * fx < 0.0) {
            b = x;
        } else {
            a = x;
            fa = fx;
        }
        N++;
    }
    return x;
}


double chordMethod(double a, double b, double eps, int &N, double delta_round = 0.0) {
    double fa = f(a);
    double fb = f(b);
    double x, fx;

    N = 0;
    if (fa == 0.0) return a;
    if (fb == 0.0) return b;

    do {
        x = a - (b - a) * fa / (fb - fa);
        fx = f(x);
        if (delta_round > 0.0) fx = roundValue(fx, delta_round);

        if (fx == 0.0) return x;
        if (fa * fx < 0.0) {
            b = x;
            fb = fx;
        } else {
            a = x;
            fa = fx;
        }
        N++;
    } while (fabs(fx) >= eps);

    return x;
}

int main() {
    cout << "========================================\n";
    cout << "  Lab Work #1\n";
    cout << "  Variant 16: f(x) = arcsin(2x/(1+x^2)) - exp(-x^2)\n";
    cout << "  Student: Rumyantsev I.V., gr. 4343\n";
    cout << "========================================\n\n";

    cout << "Select mode:\n";
    cout << "  1 - Manual input (single method)\n";
    cout << "  2 - Research (auto data collection to CSV)\n";
    cout << "Your choice: ";

    int mode;
    cin >> mode;


    if (mode == 1) {
        cout << "\nMethods:\n";
        cout << "  1 - Bisection\n";
        cout << "  2 - Chord\n";
        cout << "Your choice: ";
        int method;
        cin >> method;

        double a = 0.0, b = 1.0;
        double eps, delta;

        cout << "Enter precision eps: ";
        cin >> eps;
        cout << "Enter delta (0 = no noise): ";
        cin >> delta;

        int iter = 0;
        double root;

        if (method == 1) {
            root = bisectionMethod(a, b, eps, iter, delta);
            cout << "\n[Bisection]\n";
        } else if (method == 2) {
            root = chordMethod(a, b, eps, iter, delta);
            cout << "\n[Chord]\n";
        } else {
            cout << "Invalid method\n";
            return 1;
        }

        cout << "  Root x* = " << setprecision(10) << root << "\n";
        cout << "  f(x*)   = " << f(root) << "\n";
        cout << "  Iters   = " << iter << "\n";
    }


    else if (mode == 2) {
        double a = 0.0, b = 1.0;

        vector<double> eps_values = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6};
        vector<double> delta_values = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5};

        string methods_names[2] = {"bisection", "chord"};

        for (int m = 0; m < 2; m++) {
            string fname_eps = methods_names[m] + "_iterations_vs_eps.csv";
            string fname_delta = methods_names[m] + "_sensitivity_delta.csv";

            ofstream ofs_eps(fname_eps);
            ofstream ofs_delta(fname_delta);

            ofs_eps << "eps;iterations;root;f_root\n";
            ofs_delta << "delta;iterations;root;f_root\n";


            for (double eps : eps_values) {
                int iter = 0;
                double root = (m == 0) ? bisectionMethod(a, b, eps, iter, 0.0)
                                       : chordMethod(a, b, eps, iter, 0.0);
                double f_root = f(root);
                ofs_eps << scientific << setprecision(6) << eps << ";"
                        << iter << ";"
                        << setprecision(10) << root << ";"
                        << f_root << "\n";
            }


            double eps_fixed = 1e-6;
            for (double delta : delta_values) {
                int iter = 0;
                double root = (m == 0) ? bisectionMethod(a, b, eps_fixed, iter, delta)
                                       : chordMethod(a, b, eps_fixed, iter, delta);
                double f_root = f(root);
                ofs_delta << scientific << setprecision(6) << delta << ";"
                          << iter << ";"
                          << setprecision(10) << root << ";"
                          << f_root << "\n";
            }

            ofs_eps.close();
            ofs_delta.close();

            cout << "Files created:\n";
            cout << "  " << fname_eps << "\n";
            cout << "  " << fname_delta << "\n\n";
        }

        cout << "Research finished. Run Python script to build plots.\n";
    }

    else {
        cout << "Invalid mode\n";
        return 1;
    }

    return 0;
}