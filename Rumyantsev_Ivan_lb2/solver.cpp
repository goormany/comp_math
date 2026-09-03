#include "solver.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <algorithm>

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

bool gauss(vector<vector<double>> A, vector<double> b, vector<double>& x) {
    int n = A.size();
    x.resize(n);
    vector<vector<double>> M(n, vector<double>(n + 1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) M[i][j] = A[i][j];
        M[i][n] = b[i];
    }

    // Прямой ход с частичным выбором ведущего элемента
    for (int k = 0; k < n - 1; ++k) {
        // Поиск ведущего элемента
        int maxRow = k;
        double maxVal = fabs(M[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (fabs(M[i][k]) > maxVal) {
                maxVal = fabs(M[i][k]);
                maxRow = i;
            }
        }
        if (maxVal < 1e-15) return false;
        if (maxRow != k) swap(M[k], M[maxRow]);

        // Обнуление
        for (int i = k + 1; i < n; ++i) {
            double factor = M[i][k] / M[k][k];
            for (int j = k; j <= n; ++j)
                M[i][j] -= factor * M[k][j];
        }
    }

    if (fabs(M[n-1][n-1]) < 1e-15) return false;

    // Обратный ход
    for (int i = n - 1; i >= 0; --i) {
        double sum = M[i][n];
        for (int j = i + 1; j < n; ++j) sum -= M[i][j] * x[j];
        x[i] = sum / M[i][i];
    }
    return true;
}

void runTest(const vector<vector<double>>& A, const vector<double>& b, const string& label) {
    int n = A.size();

    cout << "  Gauss... " << flush;
    vector<double> x_gauss;
    auto start_gauss = high_resolution_clock::now();
    bool success = gauss(A, b, x_gauss);
    auto end_gauss = high_resolution_clock::now();
    auto dur_gauss = duration_cast<microseconds>(end_gauss - start_gauss).count();

    if (!success) {
        cout << "\n" << label << ": Gauss failed (singular matrix)\n";
        return;
    }
    cout << "done (" << dur_gauss / 1000.0 << " ms)\n";

    cout << "  LAPACK... " << flush;
    vector<double> A_lapack(n * n);
    vector<double> b_lapack = b;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            A_lapack[j * n + i] = A[i][j];

    vector<int> ipiv(n);
    int nrhs = 1, lda = n, ldb = n, info = 0;

    auto start_lap = high_resolution_clock::now();
    dgesv_(&n, &nrhs, A_lapack.data(), &lda, ipiv.data(), b_lapack.data(), &ldb, &info);
    auto end_lap = high_resolution_clock::now();
    auto dur_lap = duration_cast<microseconds>(end_lap - start_lap).count();

    if (info != 0) {
        cout << "\n" << label << ": LAPACK failed (info = " << info << ")\n";
        return;
    }
    cout << "done (" << dur_lap / 1000.0 << " ms)\n";

    double norm_b = norm(b);
    double rel_gauss = norm(residual(A, x_gauss, b)) / norm_b;
    double rel_lap = norm(residual(A, b_lapack, b)) / norm_b;

    double diff_gl = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = x_gauss[i] - b_lapack[i];
        diff_gl += d * d;
    }
    diff_gl = sqrt(diff_gl);

    double diff_exact = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = x_gauss[i] - (i + 1);
        diff_exact += d * d;
    }
    diff_exact = sqrt(diff_exact);

    cout << "\n========================================\n";
    cout << "  " << label << " (" << n << "x" << n << ")\n";
    cout << "========================================\n";
    cout << fixed << setprecision(2);
    cout << "Gauss time:    " << dur_gauss / 1000.0 << " ms\n";
    cout << "LAPACK time:   " << dur_lap / 1000.0 << " ms\n";
    cout << scientific << setprecision(4);
    cout << "Gauss relres:  " << rel_gauss << "\n";
    cout << "LAPACK relres: " << rel_lap << "\n";
    cout << "Gauss vs LAP:  " << diff_gl << "\n";
    cout << "Gauss vs x*:   " << diff_exact << "\n";
}