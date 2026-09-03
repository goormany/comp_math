#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>

using namespace std;

const int NMAX = 5000;

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

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "использование: " << argv[0] << " <файл_с_матрицей>\n";
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
    vector<double> x;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            in >> A[i][j];
        }
        in >> b[i];
    }
    in.close();

    cout << "размер матрицы: " << n << "x" << n << "\n";
    cout << "решаем систему\n";

    auto start_time = chrono::high_resolution_clock::now();
    
    int status = GAUSS(A, b, x, n);
    
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;

    if (status == 1) {
        cout << "ошибка: Матрица вырождена!\n";
    } else if (status == 2) {
        cout << "ошибка: n < 2\n";
    } else if (status == 3) {
        cout << "ошибка: n > NMAX (" << NMAX << ")\n";
    } else if (status == 0) {
        cout << "система успешно решена за " << elapsed.count() << " сек.\n";

        double norm_L2 = 0.0;
        double norm_Linf = 0.0;

        for (int i = 0; i < n; ++i) {
            double ax_val = 0.0;
            for (int j = 0; j < n; ++j) {
                ax_val += A[i][j] * x[j];
            }
            double residual = std::abs(ax_val - b[i]);
            
            norm_L2 += residual * residual;
            if (residual > norm_Linf) norm_Linf = residual;
        }
        norm_L2 = std::sqrt(norm_L2);

        cout << "\nРЕЗУЛЬТАТЫ\n";
        cout << scientific << setprecision(6);
        cout << "норма невязки L2: " << norm_L2 << "\n";
        cout << "норма невязки L_inf: " << norm_Linf << "\n";
        
        if (n <= 20) {
            cout << "\nвектор решений x:\n";
            for (int i = 0; i < n; ++i) {
                cout << "x[" << i << "] = " << fixed << setprecision(6) << x[i] << "\n";
            }
        } else {
            cout << "\n(матрица слишком большая для вывода вектора решений)\n";
        }
    }

    return 0;
}