#include "matrix.h"
#include <random>
#include <chrono>

vector<vector<double>> hilbertMatrix(int n) {
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = 1.0 / (i + j + 1);
    return A;
}

vector<vector<double>> randomMatrix(int n) {
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 gen(seed);
    uniform_real_distribution<double> dis(-1.0, 1.0);
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = dis(gen);
    return A;
}

vector<vector<double>> tridiagonalMatrix(int n) {
    vector<vector<double>> A(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        A[i][i] = 2.0;
        if (i > 0) A[i][i-1] = -1.0;
        if (i < n-1) A[i][i+1] = -1.0;
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