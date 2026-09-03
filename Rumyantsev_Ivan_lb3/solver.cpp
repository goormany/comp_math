#include "solver.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

extern "C" {
    void dgesv_(int* n, int* nrhs, double* A, int* lda, int* ipiv, double* B, int* ldb, int* info);
}

using namespace chrono;

double norm(const vector<double>& v) {
    double sum = 0.0;
    for (double x : v) sum += x * x;
    return sqrt(sum);
}

vector<double> residual(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b) {
    int n = A.size();
    vector<double> r(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < n; ++j) s += A[i][j] * x[j];
        r[i] = b[i] - s;
    }
    return r;
}

bool jacobiMethod(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x, int& iterations,
                  double eps, int maxIter) {
    int n = A.size();
    vector<double> x_old(n, 0.0);
    vector<double> x_new(n);

    for (int i = 0; i < n; ++i)
        if (fabs(A[i][i]) < 1e-12) return false;

    for (iterations = 1; iterations <= maxIter; ++iterations) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j)
                if (j != i) sum += A[i][j] * x_old[j];
            x_new[i] = (b[i] - sum) / A[i][i];
        }

        double diffNorm = 0.0;
        for (int i = 0; i < n; ++i)
            diffNorm = max(diffNorm, fabs(x_new[i] - x_old[i]));

        if (diffNorm < eps) {
            x = x_new;
            return true;
        }
        x_old = x_new;
    }
    x = x_new;
    return false;
}

void runTestJacobi(const vector<vector<double>>& A, const vector<double>& b, const string& label) {
    int n = A.size();

    cout << "  Jacobi... " << flush;
    vector<double> x_jacobi;
    int iters;
    auto start_j = high_resolution_clock::now();
    bool ok = jacobiMethod(A, b, x_jacobi, iters);
    auto end_j = high_resolution_clock::now();
    auto dur_j = duration_cast<microseconds>(end_j - start_j).count();

    if (!ok) {
        cout << "NOT CONVERGED\n";
        return;
    }
    cout << "done (" << dur_j / 1000.0 << " ms, " << iters << " iters)\n";

    cout << "  LAPACK... " << flush;
    vector<double> A_lapack(n * n), b_lapack = b;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            A_lapack[j * n + i] = A[i][j];

    vector<int> ipiv(n);
    int nrhs = 1, lda = n, ldb = n, info = 0;

    auto start_l = high_resolution_clock::now();
    dgesv_(&n, &nrhs, A_lapack.data(), &lda, ipiv.data(), b_lapack.data(), &ldb, &info);
    auto end_l = high_resolution_clock::now();
    auto dur_l = duration_cast<microseconds>(end_l - start_l).count();

    if (info != 0) {
        cout << "LAPACK error\n";
        return;
    }
    cout << "done (" << dur_l / 1000.0 << " ms)\n";

    double nb = norm(b);
    double rel_j = norm(residual(A, x_jacobi, b)) / nb;
    double rel_l = norm(residual(A, b_lapack, b)) / nb;

    double diff_jl = 0.0, diff_ex = 0.0;
    for (int i = 0; i < n; ++i) {
        diff_jl += pow(x_jacobi[i] - b_lapack[i], 2);
        diff_ex += pow(x_jacobi[i] - (i + 1), 2);
    }

    cout << "\n========================================\n";
    cout << "  " << label << " (" << n << "x" << n << ")\n";
    cout << "========================================\n";
    cout << fixed << setprecision(2);
    cout << "Jacobi time:   " << dur_j / 1000.0 << " ms\n";
    cout << "LAPACK time:   " << dur_l / 1000.0 << " ms\n";
    cout << "Jacobi iters:  " << iters << "\n";
    cout << scientific << setprecision(4);
    cout << "Jacobi relres: " << rel_j << "\n";
    cout << "LAPACK relres: " << rel_l << "\n";
    cout << "Jacobi vs LAP: " << sqrt(diff_jl) << "\n";
    cout << "Jacobi vs x*:  " << sqrt(diff_ex) << "\n";
}