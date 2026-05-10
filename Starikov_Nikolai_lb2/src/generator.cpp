#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

vector<vector<double>> buildHilbert(int size) {
    vector<vector<double>> mat(size, vector<double>(size));
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            mat[i][j] = 1.0 / (i + j + 1);
    return mat;
}

vector<vector<double>> buildRandom(int size) {
    unsigned seedValue = chrono::system_clock::now().time_since_epoch().count();
    mt19937 randomEngine(seedValue);
    uniform_real_distribution<double> valueDist(-1.0, 1.0);
    vector<vector<double>> mat(size, vector<double>(size));
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            mat[i][j] = valueDist(randomEngine);
    return mat;
}

vector<double> makeRightSide(const vector<vector<double>>& mat) {
    int size = mat.size();
    vector<double> rhs(size, 0.0);
    for (int i = 0; i < size; ++i) {
        double accum = 0.0;
        for (int j = 0; j < size; ++j)
            accum += mat[i][j] * (j + 1);
        rhs[i] = accum;
    }
    return rhs;
}

int main() {
    int matrixSize, choice;
    cout << "Введите размер матрицы n: ";
    cin >> matrixSize;
    cout << "Выберите тип матрицы:\n";
    cout << "1 - Гильберт\n";
    cout << "2 - Случайная\n";
    cout << "Ваш выбор: ";
    cin >> choice;

    vector<vector<double>> matrixData;
    switch (choice) {
        case 1:
            matrixData = buildHilbert(matrixSize);
            break;
        case 2:
            matrixData = buildRandom(matrixSize);
            break;
        default:
            cerr << "Неверный тип матрицы\n";
            return 1;
    }

    vector<double> rhsVector = makeRightSide(matrixData);

    ofstream outputFile("input.txt");
    outputFile << scientific << setprecision(16);
    outputFile << matrixSize << "\n";
    for (int i = 0; i < matrixSize; ++i) {
        for (int j = 0; j < matrixSize; ++j)
            outputFile << matrixData[i][j] << " ";
        outputFile << "\n";
    }
    for (int i = 0; i < matrixSize; ++i)
        outputFile << rhsVector[i] << (i < matrixSize - 1 ? " " : "\n");
    outputFile.close();

    cout << "Матрица и вектор записаны в файл input.txt\n";
    return 0;
}