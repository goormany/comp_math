#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <random>

// LAPACK
extern "C" {
    void dgesv_(int* n, int* nrhs, double* A, int* lda, int* ipiv,
                double* B, int* ldb, int* info);
}

using namespace std;
using namespace chrono;

vector<vector<double>> randomDiagonallyDominantMatrix(int n) {
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

double norm(const vector<double>& v) {
    double sum = 0.0;
    for (double x : v) sum += x * x;
    return sqrt(sum);
}

vector<double> residual(const vector<vector<double>>& A,
                        const vector<double>& x,
                        const vector<double>& b) {
    int n = A.size();
    vector<double> r(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < n; ++j)
            s += A[i][j] * x[j];
        r[i] = b[i] - s;
    }
    return r;
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

vector<double> jacobiMethod(const vector<vector<double>>& A,
                            const vector<double>& b,
                            double epsilon, int maxIter,
                            int& iterCount, bool& converged) {
    int n = A.size();
    vector<double> x_old(n, 0.0);
    vector<double> x_new(n);
    for (int i = 0; i < n; ++i) {
        if (fabs(A[i][i]) < 1e-12) {
            converged = false;
            return {};
        }
    }
    for (iterCount = 1; iterCount <= maxIter; ++iterCount) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j)
                if (j != i) sum += A[i][j] * x_old[j];
            x_new[i] = (b[i] - sum) / A[i][i];
        }
        double diffNorm = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = fabs(x_new[i] - x_old[i]);
            if (diff > diffNorm) diffNorm = diff;
        }
        if (diffNorm < epsilon) {
            converged = true;
            return x_new;
        }
        x_old = x_new;
    }
    converged = false;
    return x_new;
}

int main() {
    int n;
    cout << "Введите размер матрицы n: ";
    cin >> n;
    if (n <= 0) {
        cerr << "Размер должен быть положительным." << endl;
        return 1;
    }

    vector<vector<double>> A = randomDiagonallyDominantMatrix(n);
    vector<double> b = generateB(A);

    vector<double> x_exact(n);
    for (int j = 0; j < n; ++j)
        x_exact[j] = j + 1;

    cout << "Сгенерирована система " << n << "x" << n;
    if (checkDiagonalDominance(A))
        cout << " (строгое диагональное преобладание).\n";
    else
        cout << " (диагональное преобладание отсутствует).\n";

    const double EPS = 1e-8;
    const int MAX_ITER = 50000;
    int iterCount = 0;
    bool jacobiConverged = false;

    auto start_jacobi = high_resolution_clock::now();
    vector<double> x_jacobi = jacobiMethod(A, b, EPS, MAX_ITER, iterCount, jacobiConverged);
    auto end_jacobi = high_resolution_clock::now();
    auto duration_jacobi = duration_cast<microseconds>(end_jacobi - start_jacobi).count();

    vector<double> A_lapack(n * n);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            A_lapack[j * n + i] = A[i][j];

    vector<double> x_lapack = b;
    vector<int> ipiv(n);
    int nrhs = 1, lda = n, ldb = n, info = 0;

    auto start_lapack = high_resolution_clock::now();
    dgesv_(&n, &nrhs, A_lapack.data(), &lda, ipiv.data(),
           x_lapack.data(), &ldb, &info);
    auto end_lapack = high_resolution_clock::now();
    auto duration_lapack = duration_cast<microseconds>(end_lapack - start_lapack).count();

    if (info != 0) {
        cerr << "LAPACK dgesv завершился с ошибкой (info = " << info << ")\n";
        return 1;
    }

    cout << fixed << setprecision(6);
    cout << "\nРазмерность системы: " << n << "\n";
    cout << "Время метода Якоби: " << duration_jacobi << " мкс";
    if (!jacobiConverged) cout << " (НЕ сошёлся)";
    cout << "\n";
    cout << "Время LAPACK (dgesv): " << duration_lapack << " мкс\n";

    if (jacobiConverged)
        cout << "Метод Якоби сошёлся за " << iterCount << " итераций\n";
    else
        cout << "Метод Якоби НЕ сошёлся за " << MAX_ITER << " итераций\n";

    double norm_b = norm(b);
    if (!x_jacobi.empty()) {
        vector<double> r_jacobi = residual(A, x_jacobi, b);
        cout << "Относительная невязка Якоби: " << scientific << norm(r_jacobi) / norm_b << "\n";
    }

    vector<double> r_lapack = residual(A, x_lapack, b);
    cout << "Относительная невязка LAPACK: " << scientific << norm(r_lapack) / norm_b << "\n";

    cout << "\nТочное решение: x_j = j+1 (j = 0.." << n-1 << ")\n";
    if (!x_jacobi.empty()) {
        double err_jacobi = 0.0;
        for (int i = 0; i < n; ++i)
            err_jacobi += (x_jacobi[i] - x_exact[i]) * (x_jacobi[i] - x_exact[i]);
        err_jacobi = sqrt(err_jacobi);
        cout << "Погрешность ||x_jacobi - x_exact||_2 : " << scientific << err_jacobi << "\n";
    }

    double err_lapack = 0.0;
    for (int i = 0; i < n; ++i)
        err_lapack += (x_lapack[i] - x_exact[i]) * (x_lapack[i] - x_exact[i]);
    err_lapack = sqrt(err_lapack);
    cout << "Погрешность ||x_lapack - x_exact||_2 : " << scientific << err_lapack << "\n";

    return 0;
}