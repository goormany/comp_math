#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>
#include <chrono>

using namespace std;

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

int main() {
    int n, type;
    cout << "Введите размер матрицы n: ";
    cin >> n;
    cout << "Выберите тип матрицы:\n";
    cout << "1 - Гильберт (плохо обусловленная)\n";
    cout << "2 - Случайная (умеренно обусловленная)\n";
    cout << "3 - Трёхдиагональная (хорошо обусловленная)\n";
    cout << "Ваш выбор: ";
    cin >> type;

    vector<vector<double>> A;
    switch (type) {
        case 1: A = hilbertMatrix(n); break;
        case 2: A = randomMatrix(n); break;
        case 3: A = tridiagonalMatrix(n); break;
        default:
            cerr << "Неверный тип матрицы\n";
            return 1;
    }

    vector<double> b = generateB(A);

    ofstream fout("input.txt");
    fout << scientific << setprecision(16);
    fout << n << "\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            fout << A[i][j] << " ";
        fout << "\n";
    }
    for (int i = 0; i < n; ++i)
        fout << b[i] << (i < n-1 ? " " : "\n");
    fout.close();

    cout << "Матрица и вектор записаны в файл input.txt\n";
    return 0;
}