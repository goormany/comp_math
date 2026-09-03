#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <iomanip>
using namespace std;


void generateTestMatrix(int n, const string& filename, int seed) {
    srand(seed);
    
    // Генерируем случайное решение
    vector<double> x(n);
    for (int i = 0; i < n; i++) {
        x[i] = rand() % 20 + 1;
    }
    
    // Генерируем матрицу A
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                A[i][j] = (rand() % 90 + 10) + (rand() % 100) / 100.0;
            } else {
                A[i][j] = (rand() % 20 - 10) + (rand() % 100) / 100.0;
            }
        }
    }
    
    // Вычисляем вектор b = A * x
    vector<double> b(n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i] += A[i][j] * x[j];
        }
    }
    
    ofstream file(filename);
    file << n << endl;
    
    // Матрица A
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file << fixed << setprecision(6) << A[i][j];
            if (j < n - 1) file << " ";
        }
        file << endl;
    }
    
    // Вектор b
    for (int i = 0; i < n; i++) {
        file << fixed << setprecision(6) << b[i];
        if (i < n - 1) file << " ";
    }
    file << endl;
}

void generateHilbertMatrix(int n, const string& filename) {
    // Матрица Гильберта: H[i][j] = 1.0 / (i + j + 1)
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = 1.0 / (i + j + 1);
        }
    }
    
    vector<double> x(n, 1.0);
    
    // Вычисляем вектор b = A * x
    vector<double> b(n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i] += A[i][j] * x[j];  // b[i] = Σ A[i][j] * 1
        }
    }
    
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Ошибка: не удалось создать файл " << filename << endl;
        return;
    }
    
    file << n << endl;
    
    // Матрица A
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file << scientific << setprecision(15) << A[i][j];
            if (j < n - 1) file << " ";
        }
        file << endl;
    }
    
    // Вектор b
    for (int i = 0; i < n; i++) {
        file << scientific << setprecision(15) << b[i];
        if (i < n - 1) file << " ";
    }
    file << endl;
}

int main() {
    generateTestMatrix(500, "out_500.txt", 789);
    generateTestMatrix(750, "out_750.txt", 999);
    generateTestMatrix(1000, "out_1000.txt", 1342);
    generateHilbertMatrix(500, "out_500_g.txt");
    generateHilbertMatrix(750, "out_750_g.txt");
    generateHilbertMatrix(1000, "out_1000_g.txt");
    return 0;
}