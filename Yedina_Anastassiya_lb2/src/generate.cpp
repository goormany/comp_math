#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

void generateRandomDiagonallyDominant(int n, vector<vector<double>>& A, vector<double>& b) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(-50.0, 50.0);

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                A[i][j] = dist(gen);
                sum += std::abs(A[i][j]);
            }
        }
        A[i][i] = sum + std::abs(dist(gen)) + 1.0;
        if (dist(gen) < 0) A[i][i] = -A[i][i];

        b[i] = dist(gen);
    }
}

void generateHilbert(int n, vector<vector<double>>& A, vector<double>& b) {
    for (int i = 0; i < n; ++i) {
        b[i] = 0.0;
        for (int j = 0; j < n; ++j) {
            A[i][j] = 1.0 / (i + j + 1.0);
            b[i] += A[i][j];
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cout << "неверные арги\n";
        return 1;
    }

    int n = stoi(argv[1]);
    int type = stoi(argv[2]);
    string filename = argv[3];

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);

    if (type == 1) {
        generateRandomDiagonallyDominant(n, A, b);
        cout << "сгенерирована случайная матрица с диагональным преобладанием\n";
    } else if (type == 2) {
        generateHilbert(n, A, b);
        cout << "сгенерирована матрица Гильберта\n";
    } else {
        cout << "введенно неверное значение\n";
        return 1;
    }

    ofstream out(filename);
    if (!out.is_open()) {
        cerr << "ошибка создания файла!\n";
        return 1;
    }

    out << n << "\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            out << fixed << setprecision(10) << A[i][j] << " ";
        }
        out << b[i] << "\n";
    }

    out.close();
    cout << "данные записаны в " << filename << "\n";
    return 0;
}
