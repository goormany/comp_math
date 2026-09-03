#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

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
        cout << "гайд: " << argv[0] << " <файл_с_матрицей>\n";
        return 1;
    }

    ifstream in(argv[1]);
    if (!in.is_open()) {
        cerr << "ошибка: не удалось открыть файл.\n";
        return 1;
    }

    int n;
    in >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            in >> A[i][j];
        }
        in >> b[i];
    }
    in.close();

    cout << "размер матрицы: " << n << "x" << n << "\n";
    cout << "решение системы методом Якоби (простых итераций)\n";

    const int maxIter = 100000;
    const double tol = 1e-12;

    vector<double> x;
    int iter_count;
    double residual_L2, residual_Linf;

    auto start_time = chrono::high_resolution_clock::now();

    bool converged = jacobi(A, b, x, maxIter, tol, iter_count, residual_L2, residual_Linf);

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;

    if (converged) {
        cout << "метод сошёлся за " << iter_count << " итераций.\n";
    } else {
        cout << "метод не сошёлся за " << maxIter << " итераций.\n";
    }
    cout << "время выполнения: " << elapsed.count() << " сек.\n";

    cout << "\nРЕЗУЛЬТАТЫ\n";
    cout << scientific << setprecision(6);
    cout << "Норма невязки L2:   " << residual_L2 << "\n";
    cout << "Норма невязки L_inf: " << residual_Linf << "\n";

    // первые 10 компонент решения (матрица не слишком большая)
    if (n <= 20) {
        cout << "\nвектор решения x:\n";
        cout << fixed << setprecision(8);
        for (int i = 0; i < n; ++i) {
            cout << "x[" << i << "] = " << x[i] << "\n";
        }
    } else {
        cout << "\n(размер матрицы > 20, вывод решения опущен)\n";
        cout << "первые 5 компонент: ";
        for (int i = 0; i < min(5, n); ++i) cout << x[i] << " ";
        cout << "\n";
    }

    return 0;
}