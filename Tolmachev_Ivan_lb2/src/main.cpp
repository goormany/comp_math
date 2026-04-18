#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <algorithm>

extern "C" {
    void dgesv_(int* n, int* nrhs, double* A, int* lda, int* ipiv, double* B, int* ldb, int* info);
}

using namespace std;
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

    for (int k = 0; k < n - 1; ++k) {
        int maxRow = k;
        double maxVal = fabs(M[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (fabs(M[i][k]) > maxVal) {
                maxVal = fabs(M[i][k]);
                maxRow = i;
            }
        }
        if (maxVal == 0.0 || isnan(maxVal)) {
            return false;
        }
        if (maxRow != k) swap(M[k], M[maxRow]);

        for (int i = k + 1; i < n; ++i) {
            double factor = M[i][k] / M[k][k];
            for (int j = k; j <= n; ++j) {
                M[i][j] -= factor * M[k][j];
            }
        }
    }

    if (fabs(M[n-1][n-1]) == 0.0 || isnan(M[n-1][n-1])) return false;

    for (int i = n - 1; i >= 0; --i) {
        double sum = M[i][n];
        for (int j = i + 1; j < n; ++j) sum -= M[i][j] * x[j];
        x[i] = sum / M[i][i];
    }
    return true;
}

int main() {
    ifstream fin("input.txt");
    if (!fin) {
        cerr << "Не удалось открыть файл input.txt" << endl;
        return 1;
    }

    int n;
    fin >> n;
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            fin >> A[i][j];

    for (int i = 0; i < n; ++i)
        fin >> b[i];
    fin.close();

    vector<double> A_lapack(n * n);
    vector<double> b_lapack = b;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            A_lapack[j * n + i] = A[i][j];

    vector<double> x_gauss;
    auto start_gauss = high_resolution_clock::now();
    bool success = gauss(A, b, x_gauss);
    auto end_gauss = high_resolution_clock::now();
    auto duration_gauss = duration_cast<microseconds>(end_gauss - start_gauss).count();

    if (!success) {
        cerr << "Метод Гаусса: матрица вырождена или близка к вырожденной" << endl;
        return 1;
    }

    vector<int> ipiv(n);
    int nrhs = 1;
    int lda = n;
    int ldb = n;
    int info = 0;

    auto start_lapack = high_resolution_clock::now();
    dgesv_(&n, &nrhs, A_lapack.data(), &lda, ipiv.data(), b_lapack.data(), &ldb, &info);
    auto end_lapack = high_resolution_clock::now();
    auto duration_lapack = duration_cast<microseconds>(end_lapack - start_lapack).count();

    if (info != 0) {
        cerr << "LAPACK: ошибка решения (info = " << info << ")" << endl;
        return 1;
    }

    vector<double> r_gauss = residual(A, x_gauss, b);
    vector<double> r_lapack = residual(A, b_lapack, b);
    double norm_b = norm(b);
    double relres_gauss = norm(r_gauss) / norm_b;
    double relres_lapack = norm(r_lapack) / norm_b;

    double err = 0.0;
    for (int i = 0; i < n; ++i) {
        double diff = x_gauss[i] - b_lapack[i];
        err += diff * diff;
    }
    err = sqrt(err);

    cout << fixed << setprecision(6);
    cout << "Размерность системы: " << n << "\n";
    cout << "Время метода Гаусса: " << duration_gauss << " мкс\n";
    cout << "Время LAPACK (dgesv): " << duration_lapack << " мкс\n";
    cout << "Относительная невязка Гаусса: " << scientific << relres_gauss << "\n";
    cout << "Относительная невязка LAPACK: " << scientific << relres_lapack << "\n";
    cout << "Погрешность решения Гаусса относительно LAPACK (евклидова норма): " << scientific << err << "\n";

    ofstream fout("solution.txt");
    fout << scientific << setprecision(16);
    for (int i = 0; i < n; ++i)
        fout << x_gauss[i] << (i < n-1 ? " " : "\n");
    fout.close();
    cout << "Решение сохранено в файл solution.txt\n";

    return 0;
}