#include "matrix.h"
#include <random>
#include <chrono>
#include <cmath>

vector<vector<double>> diagonallyDominantMatrix(int n) {
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 gen(seed);
    uniform_real_distribution<double> dis(-1.0, 1.0);
    uniform_real_distribution<double> pos_dis(0.1, 1.0);

    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
        double off_diag_sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                A[i][j] = dis(gen);
                off_diag_sum += fabs(A[i][j]);
            }
        }
        A[i][i] = off_diag_sum + pos_dis(gen);
    }
    return A;
}

vector<double> generateB(const vector<vector<double>>& A) {
    int n = A.size();
    vector<double> b(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j)
            sum += A[i][j] * (j + 1);
        b[i] = sum;
    }
    return b;
}

bool checkDiagonalDominance(const vector<vector<double>>& A) {
    int n = A.size();
    for (int i = 0; i < n; ++i) {
        double diag = fabs(A[i][i]);
        double sum = 0.0;
        for (int j = 0; j < n; ++j)
            if (j != i) sum += fabs(A[i][j]);
        if (diag <= sum) return false;
    }
    return true;
}