#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>
#include <cstdlib>

using namespace std;

double randomDouble(double min, double max, mt19937& gen) {
    uniform_real_distribution<double> dist(min, max);
    return dist(gen);
}

void generateRandomSystem(int n, vector<vector<double>>& A, vector<double>& b,
                          double minVal, double maxVal, mt19937& gen) {
    A.assign(n, vector<double>(n, 0.0));
    b.assign(n, 0.0);

    uniform_real_distribution<double> offDiagDist(minVal, maxVal);

    for (int i = 0; i < n; ++i) {
        double offDiagSum = 0.0;

        for (int j = 0; j < n; ++j) {
            if (j != i) {
                A[i][j] = offDiagDist(gen);
                offDiagSum += fabs(A[i][j]);
            }
        }

        double diagExtra = randomDouble(1.0, maxVal/2.0, gen);
        double diagValue = offDiagSum + diagExtra;

        if (randomDouble(0.0, 1.0, gen) < 0.5)
            diagValue = -diagValue;

        A[i][i] = diagValue;
    }

    uniform_real_distribution<double> rhsDist(minVal, maxVal);
    for (int i = 0; i < n; ++i) {
        b[i] = rhsDist(gen);
    }
}

void printSystem(const vector<vector<double>>& A, const vector<double>& b) {
    int n = A.size();
    cout << "Сгенерированная расширенная матрица (A|b):\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << setw(10) << fixed << setprecision(4) << A[i][j] << " ";
        }
        cout << "| " << setw(10) << fixed << setprecision(4) << b[i] << "\n";
    }
}


void saveToFile(const string& filename, const vector<vector<double>>& A, const vector<double>& b) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Ошибка: не удалось открыть файл для записи " << filename << endl;
        exit(1);
    }

    int n = A.size();
    file << n << "\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file << A[i][j] << " ";
        }
        file << b[i] << "\n";
    }
    file.close();
    cout << "Система сохранена в файл: " << filename << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Использование: " << argv[0] << " <n> <выходной_файл> [min_value] [max_value]\n";
        cerr << "  n           - количество уравнений (целое положительное число)\n";
        cerr << "  выходной_файл - имя файла для сохранения системы\n";
        cerr << "  min_value   - минимальное значение для генерации (по умолчанию -10.0)\n";
        cerr << "  max_value   - максимальное значение для генерации (по умолчанию 10.0)\n";
        return 1;
    }

    int n = atoi(argv[1]);
    if (n <= 0) {
        cerr << "Ошибка: n должно быть положительным целым числом.\n";
        return 1;
    }

    string filename = argv[2];
    double minVal = -10.0, maxVal = 10.0;
    if (argc >= 4) minVal = atof(argv[3]);
    if (argc >= 5) maxVal = atof(argv[4]);

    if (minVal >= maxVal) {
        cerr << "Ошибка: min_value должно быть меньше max_value.\n";
        return 1;
    }

    random_device rd;
    mt19937 gen(rd());

    vector<vector<double>> A;
    vector<double> b;
    generateRandomSystem(n, A, b, minVal, maxVal, gen);

    printSystem(A, b);
    cout << endl;

    saveToFile(filename, A, b);

    return 0;
}