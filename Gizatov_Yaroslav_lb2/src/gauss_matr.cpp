#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>

using namespace std;

const double EPS = 1e-9;


void printMatrix(const vector<vector<double>>& A) {
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            cout << setw(10) << A[i][j] << " ";
        }
        cout << endl;
    }
}


bool gauss(vector<vector<double>>& A, vector<double>& x) {
    int n = A.size();

    for (int i = 0; i < n; ++i) {
        int pivot = i;
        double max_val = fabs(A[i][i]);
        for (int j = i + 1; j < n; ++j) {
            if (fabs(A[j][i]) > max_val) {
                max_val = fabs(A[j][i]);
                pivot = j;
            }
        }
        if (max_val < EPS) {
            cerr << "Ошибка: матрица вырождена." << endl;
            return false;
        }
        if (pivot != i) {
            swap(A[i], A[pivot]);
        }

        for (int j = i + 1; j < n; ++j) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k <= n; ++k) {
                A[j][k] -= factor * A[i][k];
            }
        }
    }

    x.resize(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        x[i] = (A[i][n] - sum) / A[i][i];
    }
    return true;
}

void computeResidual(const vector<vector<double>>& A_orig, const vector<double>& b_orig,
                     const vector<double>& x, double& norm_l2, double& norm_linf) {
    int n = A_orig.size();
    norm_l2 = 0.0;
    norm_linf = 0.0;
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += A_orig[i][j] * x[j];
        }
        double resid = sum - b_orig[i];
        norm_l2 += resid * resid;
        if (fabs(resid) > norm_linf) norm_linf = fabs(resid);
    }
    norm_l2 = sqrt(norm_l2);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Использование: " << argv[0] << " <входной_файл> <выходной_файл>" << endl;
        return 1;
    }

    string input_file = argv[1];
    string output_file = argv[2];

    ifstream fin(input_file);
    if (!fin.is_open()) {
        cerr << "Ошибка: не удалось открыть входной файл " << input_file << endl;
        return 1;
    }

    int n;
    fin >> n;
    if (fin.fail() || n <= 0) {
        cerr << "Ошибка: первая строка должна содержать положительное целое число." << endl;
        return 1;
    }

    vector<vector<double>> A(n, vector<double>(n + 1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (!(fin >> A[i][j])) {
                cerr << "Ошибка: недостаточно данных." << endl;
                return 1;
            }
        }
    }
    fin.close();

    vector<vector<double>> A_orig(n, vector<double>(n));
    vector<double> b_orig(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) A_orig[i][j] = A[i][j];
        b_orig[i] = A[i][n];
    }

    auto start = chrono::high_resolution_clock::now();
    vector<double> x;
    bool success = gauss(A, x);
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double>(end - start).count();

    if (!success) {
        cerr << "Решение не найдено (матрица вырождена)." << endl;
        return 1;
    }

    double norm_l2, norm_linf;
    computeResidual(A_orig, b_orig, x, norm_l2, norm_linf);

    ofstream fout(output_file);
    if (!fout.is_open()) {
        cerr << "Ошибка: не удалось создать выходной файл " << output_file << endl;
        return 1;
    }

    fout << fixed << setprecision(12);
    for (int i = 0; i < n; ++i) {
        fout << "x" << i + 1 << " = " << x[i] << endl;
    }
    fout << "\n# Невязка (L2) = " << norm_l2 << endl;
    fout << "# Невязка (L∞) = " << norm_linf << endl;
    fout << "# Время решения = " << elapsed << " сек." << endl;
    fout.close();

    cout << "Решение записано в файл: " << output_file << endl;
    cout << "Размер системы: " << n << " x " << n << endl;
    cout << "Невязка L2 = " << norm_l2 << ", L∞ = " << norm_linf << endl;
    cout << "Время решения: " << elapsed << " сек." << endl;

    return 0;
}