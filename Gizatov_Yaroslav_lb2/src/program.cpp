#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <random>
#include <windows.h>
#include <clocale>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

long double EPS = 1e-12L;


bool gauss_complete_pivoting(vector<vector<long double>>& A, vector<long double>& x) {
    int n = A.size();
    vector<int> row_perm(n), col_perm(n);
    for (int i = 0; i < n; ++i) row_perm[i] = col_perm[i] = i;

    for (int k = 0; k < n; ++k) {
        int max_row = k, max_col = k;
        long double max_val = fabsl(A[k][k]);
        for (int i = k; i < n; ++i)
            for (int j = k; j < n; ++j)
                if (fabsl(A[i][j]) > max_val) {
                    max_val = fabsl(A[i][j]);
                    max_row = i; max_col = j;
                }
        if (max_val < 1e-14L) {
            cerr << "Ошибка: матрица вырождена (Гаусс).\n";
            return false;
        }
        if (max_row != k) {
            swap(A[k], A[max_row]);
            swap(row_perm[k], row_perm[max_row]);
        }
        if (max_col != k) {
            for (int i = 0; i < n; ++i) swap(A[i][k], A[i][max_col]);
            swap(col_perm[k], col_perm[max_col]);
        }
        for (int i = k + 1; i < n; ++i) {
            long double factor = A[i][k] / A[k][k];
            for (int j = k; j <= n; ++j)
                A[i][j] -= factor * A[k][j];
        }
    }

    vector<long double> y(n);
    for (int i = n - 1; i >= 0; --i) {
        long double sum = 0.0L;
        for (int j = i + 1; j < n; ++j) sum += A[i][j] * y[j];
        y[i] = (A[i][n] - sum) / A[i][i];
    }
    x.resize(n);
    for (int i = 0; i < n; ++i) x[col_perm[i]] = y[i];
    return true;
}


void computeResidual(const vector<vector<long double>>& A_orig,
                     const vector<long double>& b_orig,
                     const vector<long double>& x,
                     long double& norm_l2, long double& norm_linf) {
    int n = A_orig.size();
    norm_l2 = 0.0L;
    norm_linf = 0.0L;
    for (int i = 0; i < n; ++i) {
        long double sum = 0.0L;
        for (int j = 0; j < n; ++j) sum += A_orig[i][j] * x[j];
        long double resid = sum - b_orig[i];
        norm_l2 += resid * resid;
        if (fabsl(resid) > norm_linf) norm_linf = fabsl(resid);
    }
    norm_l2 = sqrtl(norm_l2);
}


void generateRandomSystem(int n, MatrixXd& A, VectorXd& b, VectorXd& x_exact) {
    A = MatrixXd::Random(n, n).array() * 10.0;
    x_exact = VectorXd::Random(n).array() * 10.0;
    b = A * x_exact;
}


void generateDiagonallyDominantSystem(int n, MatrixXd& A, VectorXd& b, VectorXd& x_exact) {
    A = MatrixXd::Random(n, n).array() * 10.0;

    for (int i = 0; i < n; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < n; ++j) if (j != i) row_sum += fabs(A(i, j));
        A(i, i) = row_sum + fabs(A(i, i)) + 1.0;
    }
    x_exact = VectorXd::Random(n).array() * 10.0;
    b = A * x_exact;
}


void generateHilbertSystem(int n, MatrixXd& A, VectorXd& b, VectorXd& x_exact) {
    A = MatrixXd::Zero(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A(i, j) = 1.0 / (i + j + 1);
    x_exact = VectorXd::Random(n).array() * 10.0;
    b = A * x_exact;
}


void generateSingularSystem(int n, MatrixXd& A, VectorXd& b) {
    A = MatrixXd::Random(n, n).array() * 10.0;
    if (n >= 3) {
        A.row(n-1) = A.row(0) + A.row(1);
    } else if (n == 2) {
        A.row(1) = A.row(0) * 2.0;
    } else if (n == 1) {
        A(0,0) = 0.0;
    }
    b = VectorXd::Random(n).array() * 10.0;
    FullPivLU<MatrixXd> lu(A);
    cout << "Ранг матрицы: " << lu.rank() << " из " << n << " (вырождена)\n";
    cout << "Определитель: " << scientific << A.determinant() << " (близок к 0)\n";
}


bool solveEigen(const MatrixXd& A, const VectorXd& b, VectorXd& x, double& elapsed, double& l2, double& linf) {
    auto start = chrono::high_resolution_clock::now();
    FullPivLU<MatrixXd> lu(A);
    bool invertible = lu.isInvertible();
    if (!invertible) {
        elapsed = 0.0;
        l2 = linf = 0.0;
        return false;
    }
    x = lu.solve(b);
    auto end = chrono::high_resolution_clock::now();
    elapsed = chrono::duration<double>(end - start).count();
    VectorXd residual = A * x - b;
    l2 = residual.norm();
    linf = residual.lpNorm<Infinity>();
    return true;
}


bool solveGauss(const MatrixXd& A_eigen, const VectorXd& b_eigen, vector<long double>& x,
                double& elapsed, long double& l2, long double& linf) {
    int n = A_eigen.rows();
    vector<vector<long double>> A(n, vector<long double>(n + 1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) A[i][j] = (long double)A_eigen(i, j);
        A[i][n] = (long double)b_eigen(i);
    }
    vector<vector<long double>> A_orig(n, vector<long double>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) A_orig[i][j] = A[i][j];
    vector<long double> b_orig(n);
    for (int i = 0; i < n; ++i) b_orig[i] = (long double)b_eigen(i);

    auto start = chrono::high_resolution_clock::now();
    bool success = gauss_complete_pivoting(A, x);
    auto end = chrono::high_resolution_clock::now();
    elapsed = chrono::duration<double>(end - start).count();
    if (success) computeResidual(A_orig, b_orig, x, l2, linf);
    return success;
}


double computeJacobiNorm(const MatrixXd& A) {
    int n = A.rows();
    double max_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        double diag = fabs(A(i,i));
        if (diag < 1e-15) return 1e20;
        for (int j = 0; j < n; ++j) {
            if (j != i) sum += fabs(A(i,j)) / diag;
        }
        if (sum > max_sum) max_sum = sum;
    }
    return max_sum;
}

bool jacobi(const MatrixXd& A, const VectorXd& b, VectorXd& x,
            double& elapsed, double& l2, double& linf,
            int maxIter = 10000, double tol = EPS) {
    tol = EPS;
    int n = A.rows();
    for (int i = 0; i < n; ++i) {
        if (fabs(A(i,i)) < 1e-15) {
            cerr << "Метод Якоби: нулевой диагональный элемент, метод неприменим.\n";
            return false;
        }
    }

    double normB = computeJacobiNorm(A);
    cout << "Норма матрицы Якоби ||B||_∞ = " << normB << "\n";
    if (normB >= 1.0) {
        cout << "Предупреждение: ||B|| >= 1, метод может расходиться или сходиться медленно.\n";
    }

    VectorXd x_new = VectorXd::Zero(n);
    x = VectorXd::Zero(n);
    auto start = chrono::high_resolution_clock::now();

    for (int iter = 0; iter < maxIter; ++iter) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) sum += A(i, j) * x(j);
            }
            x_new(i) = (b(i) - sum) / A(i, i);
        }
        double err = (x_new - x).norm();
        x = x_new;
        if (err < tol) {
            auto end = chrono::high_resolution_clock::now();
            elapsed = chrono::duration<double>(end - start).count();
            VectorXd residual = A * x - b;
            l2 = residual.norm();
            linf = residual.lpNorm<Infinity>();
            cout << "Метод Якоби сошёлся за " << iter+1 << " итераций.\n";
            return true;
        }
    }
    auto end = chrono::high_resolution_clock::now();
    elapsed = chrono::duration<double>(end - start).count();
    VectorXd residual = A * x - b;
    l2 = residual.norm();
    linf = residual.lpNorm<Infinity>();
    cout << "Метод Якоби не сошёлся за " << maxIter << " итераций.\n";
    return false;
}


int main() {
    SetConsoleOutputCP(CP_UTF8);
    setlocale(LC_ALL, ".UTF-8");

    MatrixXd A;
    VectorXd b, x_exact;
    int n = 0;
    bool systemReady = false;
    bool hasExact = false;

    int choice;
    do {
        cout << "\n         МЕНЮ          \n";
        cout << "1. Решить библиотечным методом (Eigen)\n";
        cout << "2. Решить методом Гаусса (полный выбор)\n";
        cout << "3. Решить методом Якоби (простых итераций)\n";
        cout << "4. Задать eps для метода Якоби\n\n";
    
        cout << "5. Сравнить решения (разность векторов и эталон)\n\n";
        
        cout << "6. Сгенерировать случайную СЛАУ (n x n)\n";        
        cout << "7. Сгенерировать матрицу Гильберта (n x n)\n";
        cout << "8. Сгенерировать ВЫРОЖДЕННУЮ матрицу (тест реакции методов)\n";
        cout << "9. Сгенерировать диагонально преобладающую СЛАУ (для Якоби)\n";
        cout << "0. Выход\n";
        cout << "Ваш выбор: ";

        cin >> choice;

        switch (choice) {
            case 6: {
                cout << "Введите размер матрицы n: ";
                cin >> n;
                if (n <= 0) {
                    cerr << "Неверный размер.\n";
                    break;
                }
                generateRandomSystem(n, A, b, x_exact);
                systemReady = true;
                hasExact = true;
                cout << "Сгенерирована случайная система " << n << "x" << n << ".\n";
                cout << "Пример: A(0,0) = " << A(0,0) << ", b(0) = " << b(0) << endl;
                break;
            }
            case 1: {
                if (!systemReady) {
                    cout << "Сначала сгенерируйте систему (пункт 1,5,6 или 8).\n";
                    break;
                }
                VectorXd x_eigen;
                double elapsed, l2, linf;
                bool ok = solveEigen(A, b, x_eigen, elapsed, l2, linf);
                cout << fixed << setprecision(12);
                cout << "\n Решение Eigen \n";
                if (ok) {
                    cout << "Матрица невырождена.\n";
                    cout << "Невязка L2   = " << scientific << l2 << "\n";
                    cout << "Невязка L∞   = " << linf << "\n";
                    cout << "Время решения = " << elapsed << " сек.\n";
                    if (hasExact) {
                        VectorXd err = x_eigen - x_exact;
                        cout << "|| x - x* ||_2 = " << err.norm() << "\n";
                        cout << "|| x - x* ||_∞ = " << err.lpNorm<Infinity>() << "\n";
                    }
                    static VectorXd x_eigen_global;
                    x_eigen_global = x_eigen;
                    cout << "Первые 5 компонент решения:\n";
                    for (int i = 0; i < min(5, n); ++i)
                        cout << "x[" << i+1 << "] = " << x_eigen(i) << "\n";
                } else {
                    cout << "Матрица вырождена, решение невозможно.\n";
                }
                break;
            }
            case 2: {
                if (!systemReady) {
                    cout << "Сначала сгенерируйте систему (пункт 1,5,6 или 8).\n";
                    break;
                }
                vector<long double> x_gauss;
                double elapsed;
                long double l2, linf;
                bool ok = solveGauss(A, b, x_gauss, elapsed, l2, linf);
                cout << fixed << setprecision(12);
                cout << "\n Решение Гаусса (полный выбор) \n";
                if (ok) {
                    cout << "Невязка L2   = " << scientific << (double)l2 << "\n";
                    cout << "Невязка L∞   = " << (double)linf << "\n";
                    cout << "Время решения = " << elapsed << " сек.\n";
                    if (hasExact) {
                        VectorXd x_gauss_eigen(n);
                        for (int i = 0; i < n; ++i) x_gauss_eigen(i) = (double)x_gauss[i];
                        VectorXd err = x_gauss_eigen - x_exact;
                        cout << "|| x - x* ||_2 = " << err.norm() << "\n";
                        cout << "|| x - x* ||_∞ = " << err.lpNorm<Infinity>() << "\n";
                    }
                    static vector<long double> x_gauss_global;
                    x_gauss_global = x_gauss;
                    cout << "Первые 5 компонент решения:\n";
                    for (int i = 0; i < min(5, n); ++i)
                        cout << "x[" << i+1 << "] = " << (double)x_gauss[i] << "\n";
                } else {
                    cout << "Метод Гаусса не смог решить (матрица вырождена или ошибка округления).\n";
                }
                break;
            }
            case 5: {
                if (!systemReady) {
                    cout << "Сначала сгенерируйте систему (пункт 1,5,6 или 8).\n";
                    break;
                }
                cout << "Хотите выполнить сравнение? (y/n): ";
                char ans; cin >> ans;
                if (ans != 'y' && ans != 'Y') break;

                VectorXd x_eigen;
                double el1, l2e, linfe;
                bool eigen_ok = solveEigen(A, b, x_eigen, el1, l2e, linfe);
                vector<long double> x_gauss_vec;
                double el2;
                long double l2g, linfg;
                bool gauss_ok = solveGauss(A, b, x_gauss_vec, el2, l2g, linfg);
                VectorXd x_jacobi;
                double el3, l2j, linfj;
                bool jacobi_ok = jacobi(A, b, x_jacobi, el3, l2j, linfj);

                cout << "\nСРАВНЕНИЕ РЕШЕНИЙ \n";
                cout << scientific << setprecision(12);

                if (hasExact) {
                    cout << "\nПогрешность относительно точного решения\n";
                    if (eigen_ok) {
                        VectorXd err_e = x_eigen - x_exact;
                        cout << "Eigen:  ||x-x*||_2 = " << err_e.norm() 
                            << ",   ||x-x*||_∞ = " << err_e.lpNorm<Infinity>() << "\n";
                    }
                    if (gauss_ok) {
                        VectorXd x_gauss_eigen(n);
                        for (int i = 0; i < n; ++i) x_gauss_eigen(i) = (double)x_gauss_vec[i];
                        VectorXd err_g = x_gauss_eigen - x_exact;
                        cout << "Gauss:  ||x-x*||_2 = " << err_g.norm() 
                            << ",   ||x-x*||_∞ = " << err_g.lpNorm<Infinity>() << "\n";
                    }
                    if (jacobi_ok) {
                        VectorXd err_j = x_jacobi - x_exact;
                        cout << "Jacobi: ||x-x*||_2 = " << err_j.norm() 
                            << ",   ||x-x*||_∞ = " << err_j.lpNorm<Infinity>() << "\n";
                    }
                }

                cout << "\nРазности между методами\n";
                if (eigen_ok && gauss_ok) {
                    VectorXd x_gauss_eigen(n);
                    for (int i = 0; i < n; ++i) x_gauss_eigen(i) = (double)x_gauss_vec[i];
                    VectorXd diff_eg = x_eigen - x_gauss_eigen;
                    cout << "||Eigen - Gauss||_2 = " << diff_eg.norm() 
                        << ",   ||Eigen - Gauss||_∞ = " << diff_eg.lpNorm<Infinity>() << "\n";
                }
                if (eigen_ok && jacobi_ok) {
                    VectorXd diff_ej = x_eigen - x_jacobi;
                    cout << "||Eigen - Jacobi||_2 = " << diff_ej.norm() 
                        << ",   ||Eigen - Jacobi||_∞ = " << diff_ej.lpNorm<Infinity>() << "\n";
                }
                if (gauss_ok && jacobi_ok) {
                    VectorXd x_gauss_eigen(n);
                    for (int i = 0; i < n; ++i) x_gauss_eigen(i) = (double)x_gauss_vec[i];
                    VectorXd diff_gj = x_gauss_eigen - x_jacobi;
                    cout << "||Gauss - Jacobi||_2 = " << diff_gj.norm() 
                        << ",   ||Gauss - Jacobi||_∞ = " << diff_gj.lpNorm<Infinity>() << "\n";
                }
                break;
            }
            case 7: {
                cout << "Введите размер матрицы Гильберта n: ";
                cin >> n;
                if (n <= 0) {
                    cerr << "Неверный размер.\n";
                    break;
                }
                generateHilbertSystem(n, A, b, x_exact);
                systemReady = true;
                hasExact = true;
                cout << "Сгенерирована матрица Гильберта " << n << "x" << n << " (очень плохо обусловлена).\n";
                cout << "Пример: A(0,0) = " << A(0,0) << ", b(0) = " << b(0) << endl;
                break;
            }
            case 8: {
                cout << "Введите размер матрицы n: ";
                cin >> n;
                if (n <= 0) {
                    cerr << "Неверный размер.\n";
                    break;
                }
                generateSingularSystem(n, A, b);
                systemReady = true;
                hasExact = false;
                cout << "Сгенерирована ВЫРОЖДЕННАЯ система " << n << "x" << n << ".\n";
                break;
            }
            case 3: {
                if (!systemReady) {
                    cout << "Сначала сгенерируйте систему (пункт 1,5,6 или 8).\n";
                    break;
                }
                VectorXd x_jacobi;
                double elapsed, l2, linf;
                bool ok = jacobi(A, b, x_jacobi, elapsed, l2, linf);
                cout << fixed << setprecision(12);
                cout << "\nРешение методом Якоби \n";
                if (ok) {
                    cout << "Невязка L2   = " << scientific << l2 << "\n";
                    cout << "Невязка L∞   = " << linf << "\n";
                    cout << "Время решения = " << elapsed << " сек.\n";
                    if (hasExact) {
                        VectorXd err = x_jacobi - x_exact;
                        cout << "|| x - x* ||_2 = " << err.norm() << "\n";
                        cout << "|| x - x* ||_∞ = " << err.lpNorm<Infinity>() << "\n";
                    }
                    cout << "Первые 5 компонент решения:\n";
                    for (int i = 0; i < min(5, n); ++i)
                        cout << "x[" << i+1 << "] = " << x_jacobi(i) << "\n";
                } else {
                    cout << "Метод Якоби не применим или не сошёлся.\n";
                }
                break;
            }
            case 9: {
                cout << "Введите размер матрицы n: ";
                cin >> n;
                if (n <= 0) {
                    cerr << "Неверный размер.\n";
                    break;
                }
                generateDiagonallyDominantSystem(n, A, b, x_exact);
                systemReady = true;
                hasExact = true;
                cout << "Сгенерирована диагонально преобладающая система " << n << "x" << n << ".\n";
                cout << "Пример: A(0,0) = " << A(0,0) << ", b(0) = " << b(0) << endl;
                break;
            }
            case 4: {
                cout << "Введите новый eps (Пример 1e-12)" << endl;
                cin >> EPS;
                break;
            }
            case 0:
                cout << "Выход.\n";
                break;
            default:
                cout << "Неверный пункт.\n";
        }
    } while (choice != 0);

    return 0;
}