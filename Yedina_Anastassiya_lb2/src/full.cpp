#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;


const int NMAX = 5000;


void computeResidual(const vector<vector<double>>& A, const vector<double>& b,
                     const vector<double>& x, double& l2, double& linf) {
    int n = A.size();
    l2 = 0.0;
    linf = 0.0;
    for (int i = 0; i < n; ++i) {
        double Ax = 0.0;
        for (int j = 0; j < n; ++j) Ax += A[i][j] * x[j];
        double resid = Ax - b[i];
        l2 += resid * resid;
        if (fabs(resid) > linf) linf = fabs(resid);
    }
    l2 = sqrt(l2);
}

void computeDiff(const vector<double>& x, const vector<double>& y,
                 double& diff_l2, double& diff_linf) {
    int n = x.size();
    diff_l2 = 0.0;
    diff_linf = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = x[i] - y[i];
        diff_l2 += d * d;
        if (fabs(d) > diff_linf) diff_linf = fabs(d);
    }
    diff_l2 = sqrt(diff_l2);
}


int GAUSS(const vector<vector<double>>& a_in, const vector<double>& b_in, vector<double>& x, int n) {
    if (n < 2) return 2;
    if (n > NMAX) return 3;


    vector<vector<double>> a = a_in;
    vector<double> b = b_in;
    x.assign(n, 0.0);

    for (int i = 0; i < n; ++i) {

        int pivot_row = i;
        double max_val = std::abs(a[i][i]);
        
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(a[k][i]) > max_val) {
                max_val = std::abs(a[k][i]);
                pivot_row = k;
            }
        }

        if (max_val < 1e-12) {
            return 1; 
        }

        if (pivot_row != i) {
            swap(a[i], a[pivot_row]);
            swap(b[i], b[pivot_row]);
        }

        for (int k = i + 1; k < n; ++k) {
            double factor = a[k][i] / a[i][i];
            for (int j = i; j < n; ++j) {
                a[k][j] -= factor * a[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }

    return 0;
}

// норма матрицы итераций Якоби ||B||_бесконечность
// B_ij = -a_ij / a_ii для i != j, и 0 на диагонали
double computeJacobiNorm(const vector<vector<double>>& A) {
    int n = A.size();
    double max_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double diag = fabs(A[i][i]);
        if (diag < 1e-15) return 1e20; // метод не применим
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (j != i) sum += fabs(A[i][j]) / diag;
        }
        if (sum > max_sum) max_sum = sum;
    }
    return max_sum;
}

bool jacobi(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x, int maxIter, double tol,
            int& iter_count, double& residual_norm_L2, double& residual_norm_Linf) {
    
    int n = A.size();


    for (int i = 0; i < n; ++i) {
        if (fabs(A[i][i]) < 1e-15) {
            cerr << "ошибка: нулевой диагональный элемент, метод Якоби неприменим\n";
            return false;
        }
    }

    
    double normB = computeJacobiNorm(A);
    cout << "норма матрицы итераций Якоби = " << normB << "\n";
    if (normB >= 1.0) {
        cout << "предупреждение: ||B|| >= 1, метод может расходиться или сходиться очень медленно\n";
    }


    x.assign(n, 0.0);
    vector<double> x_new(n, 0.0);

    for (iter_count = 0; iter_count < maxIter; ++iter_count) {

        // итерация метода
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) sum += A[i][j] * x[j]; // часть формулы, считаем сумму всех элементов матрицы умноженных на вектор приближение
            }
            x_new[i] = (b[i] - sum) / A[i][i]; // новое приблежиние
        }

        // считаем евкидову норму
        // суть -- с каждой итерацией мы находим все большее приблежение
        // соответсвенно разница между приближениями должна уменьшаться
        // при достижение достаточно маленькой разности можно считать что дальше итерировать нет смыысла решение найдено
        double err = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = x_new[i] - x[i];
            err += diff * diff;
        }
        err = sqrt(err); // ( (сумма всех x_new[i] - x[i])^2 )^1/2

        x = x_new;

        // если нормал2 меньше заданной погрешности мы заканчиваем итерации и возвращаем нынешенее приблежиние как решение
        if (err < tol) {
            // невязки
            residual_norm_L2 = 0.0;
            residual_norm_Linf = 0.0;
            for (int i = 0; i < n; ++i) {
                double Ax = 0.0;
                for (int j = 0; j < n; ++j) Ax += A[i][j] * x[j]; // скалярное произведение iй строки матрицы A на вектор x
                double resid = Ax - b[i]; // iтая компонента невязки
                residual_norm_L2 += resid * resid;
                if (fabs(resid) > residual_norm_Linf) residual_norm_Linf = fabs(resid);
            }
            residual_norm_L2 = sqrt(residual_norm_L2);
            return true;
        }
    }

    // если дошли сюда -- не сошлись за максчисло итераций итераций
    // посчитаем невязку для текущего x
    residual_norm_L2 = 0.0;
    residual_norm_Linf = 0.0;
    for (int i = 0; i < n; ++i) {
        double Ax = 0.0;
        for (int j = 0; j < n; ++j) Ax += A[i][j] * x[j];
        double resid = Ax - b[i];
        residual_norm_L2 += resid * resid;
        if (fabs(resid) > residual_norm_Linf) residual_norm_Linf = fabs(resid);
    }
    residual_norm_L2 = sqrt(residual_norm_L2);
    return false;
}


int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "использование: " << argv[0] << " <файл_с_матрицей>\n";
        return 1;
    }

    ifstream in(argv[1]);
    if (!in.is_open()) {
        cerr << "ошибка: не удалось открыть файл " << argv[1] << endl;
        return 1;
    }

    int n;
    in >> n;
    if (n <= 0) {
        cerr << "ошибка: неверный размер матрицы." << endl;
        return 1;
    }

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) in >> A[i][j];
        in >> b[i];
    }
    in.close();

    cout << fixed << setprecision(12);
    cout << "размер системы: " << n << "x" << n << "\n\n";

 
    // 1. метод Гаусс
    vector<double> x_gauss;
    auto start = chrono::high_resolution_clock::now();
    int gauss_status = GAUSS(A, b, x_gauss, n);
    auto end = chrono::high_resolution_clock::now();
    double time_gauss = chrono::duration<double>(end - start).count();

    double l2_gauss, linf_gauss;
    if (gauss_status == 0) {
        computeResidual(A, b, x_gauss, l2_gauss, linf_gauss);
    }

 
    // 2. метод Якоби
    vector<double> x_jacobi;
    int iter_count;
    double l2_jacobi, linf_jacobi;
    bool jacobi_ok = false;
    double time_jacobi = 0.0;
    if (gauss_status != 1) {
        start = chrono::high_resolution_clock::now();
        jacobi_ok = jacobi(A, b, x_jacobi, 100000, 1e-12, iter_count, l2_jacobi, linf_jacobi);
        end = chrono::high_resolution_clock::now();
        time_jacobi = chrono::duration<double>(end - start).count();
    }

 
    // 3. метод Eigen (PartialPivLU)
 
    MatrixXd A_eig(n, n);
    VectorXd b_eig(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) A_eig(i, j) = A[i][j];
        b_eig(i) = b[i];
    }
    start = chrono::high_resolution_clock::now();
    VectorXd x_eig = A_eig.partialPivLu().solve(b_eig);
    end = chrono::high_resolution_clock::now();
    double time_eigen = chrono::duration<double>(end - start).count();
    VectorXd resid_eig = A_eig * x_eig - b_eig;
    double l2_eigen = resid_eig.norm();
    double linf_eigen = resid_eig.lpNorm<Infinity>();
    vector<double> x_eigen_vec(x_eig.data(), x_eig.data() + n);

 
    cout << "метод ГАУССА\n";
    if (gauss_status == 0) {
        cout << "Время: " << time_gauss << " сек\n";
        cout << "Невязка L2   : " << scientific << l2_gauss << "\n";
        cout << "Невязка L_inf   : " << linf_gauss << "\n";
    } else if (gauss_status == 1) {
        cout << "ошибка: матрица вырождена (метод Гаусса неприменим)\n";
    } else {
        cout << "ошибка: неверный размер или превышен лимит (n<2 или n>NMAX)\n";
    }

    cout << "\nметод ЯКОБИ\n";
    if (jacobi_ok) {
        cout << "Сошёлся за " << iter_count << " итераций\n";
        cout << "Время: " << time_jacobi << " сек\n";
        cout << "Невязка L2   : " << scientific << l2_jacobi << "\n";
        cout << "Невязка L_inf   : " << linf_jacobi << "\n";
    } else {
        cout << "метод Якоби не сошёлся или неприменим (диагональные элементы близки к нулю)\n";
        if (iter_count > 0) cout << "последняя невязка L2: " << l2_jacobi << "\n";
    }

    cout << "\nметод EIGEN\n";
    cout << "Время: " << time_eigen << " сек\n";
    cout << "Невязка L2   : " << scientific << l2_eigen << "\n";
    cout << "Невязка L_inf   : " << linf_eigen << "\n";

 
    cout << "\nПОПАРНЫЕ РАЗНОСТИ РЕШЕНИЙ\n";

    if (gauss_status == 0 && jacobi_ok) {
        double diff_l2, diff_linf;
        computeDiff(x_gauss, x_jacobi, diff_l2, diff_linf);
        cout << "||Gauss - Jacobi||_2  = " << diff_l2 << "\n";
        cout << "||Gauss - Jacobi||_inf  = " << diff_linf << "\n";
    }
    if (gauss_status == 0) {
        double diff_l2, diff_linf;
        computeDiff(x_gauss, x_eigen_vec, diff_l2, diff_linf);
        cout << "||Gauss - Eigen||_2   = " << diff_l2 << "\n";
        cout << "||Gauss - Eigen||_inf   = " << diff_linf << "\n";
    }
    if (jacobi_ok) {
        double diff_l2, diff_linf;
        computeDiff(x_jacobi, x_eigen_vec, diff_l2, diff_linf);
        cout << "||Jacobi - Eigen||_2  = " << diff_l2 << "\n";
        cout << "||Jacobi - Eigen||_inf  = " << diff_linf << "\n";
    }

    if (gauss_status == 1 && !jacobi_ok) {
        cout << "матрица, вероятно, вырождена или плохо обусловлена.\n";
        cout << "результаты сравнения могут быть ненадёжными.\n";
    }

    return 0;
}