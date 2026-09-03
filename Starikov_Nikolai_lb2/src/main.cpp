#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

extern "C" {
    void dgesv_(int* n, int* nrhs, double* A, int* lda, int* ipiv, double* B, int* ldb, int* info);
}

using namespace std;
using namespace chrono;

double euclideanNorm(const vector<double>& vec) {
    double sumSq = 0.0;
    for (double val : vec) sumSq += val * val;
    return sqrt(sumSq);
}

vector<double> computeResidual(const vector<vector<double>>& mat, const vector<double>& sol, const vector<double>& rhs) {
    int size = mat.size();
    vector<double> res(size, 0.0);
    for (int i = 0; i < size; ++i) {
        double sum = 0.0;
        for (int j = 0; j < size; ++j) sum += mat[i][j] * sol[j];
        res[i] = rhs[i] - sum;
    }
    return res;
}

bool gaussElimination(vector<vector<double>> mat, vector<double> rhs, vector<double>& solution, double& minPivot, double* determinant = nullptr) {
    int size = mat.size();
    solution.resize(size);
    vector<vector<double>> extended(size, vector<double>(size + 1));
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) extended[i][j] = mat[i][j];
        extended[i][size] = rhs[i];
    }

    int swapCount = 0;
    minPivot = 1e300;

    for (int k = 0; k < size - 1; ++k) {
        int pivotRow = k;
        double maxAbs = fabs(extended[k][k]);
        for (int i = k + 1; i < size; ++i) {
            if (fabs(extended[i][k]) > maxAbs) {
                maxAbs = fabs(extended[i][k]);
                pivotRow = i;
            }
        }
        if (maxAbs == 0.0 || isnan(maxAbs)) return false;
        if (maxAbs < minPivot) minPivot = maxAbs;

        if (pivotRow != k) {
            swap(extended[k], extended[pivotRow]);
            ++swapCount;
        }

        for (int i = k + 1; i < size; ++i) {
            double factor = extended[i][k] / extended[k][k];
            for (int j = k; j <= size; ++j) {
                extended[i][j] -= factor * extended[k][j];
            }
        }
    }

    double lastPivotVal = fabs(extended[size - 1][size - 1]);
    if (lastPivotVal == 0.0 || isnan(lastPivotVal)) return false;
    if (lastPivotVal < minPivot) minPivot = lastPivotVal;

    if (determinant != nullptr) {
        *determinant = 1.0;
        for (int i = 0; i < size; ++i) *determinant *= extended[i][i];
        if (swapCount % 2 == 1) *determinant = -*determinant;
    }

    for (int i = size - 1; i >= 0; --i) {
        double sum = extended[i][size];
        for (int j = i + 1; j < size; ++j) sum -= extended[i][j] * solution[j];
        solution[i] = sum / extended[i][i];
    }
    return true;
}

int main() {
    ifstream inputFile("input.txt");
    if (!inputFile) {
        cerr << "Не удалось открыть файл input.txt" << endl;
        return 1;
    }

    int dimension;
    inputFile >> dimension;
    vector<vector<double>> matrix(dimension, vector<double>(dimension));
    vector<double> rhs(dimension);

    for (int i = 0; i < dimension; ++i)
        for (int j = 0; j < dimension; ++j)
            inputFile >> matrix[i][j];

    for (int i = 0; i < dimension; ++i)
        inputFile >> rhs[i];
    inputFile.close();

    vector<double> lapackMat(dimension * dimension);
    vector<double> lapackRhs = rhs;
    for (int j = 0; j < dimension; ++j)
        for (int i = 0; i < dimension; ++i)
            lapackMat[j * dimension + i] = matrix[i][j];

    vector<double> gaussSol;
    double smallestPivot;
    double determinant;

    auto startGauss = high_resolution_clock::now();
    bool ok = gaussElimination(matrix, rhs, gaussSol, smallestPivot, &determinant);
    auto endGauss = high_resolution_clock::now();
    auto gaussTime = duration_cast<microseconds>(endGauss - startGauss).count();

    if (!ok) {
        cerr << "Метод Гаусса: матрица вырождена (ведущий элемент = 0)" << endl;
        return 1;
    }

    const double EPS_NEAR_SINGULAR = 1e-12;
    if (smallestPivot < EPS_NEAR_SINGULAR) {
        cerr << "Предупреждение: матрица близка к вырожденной!" << endl;
        cerr << "Минимальный ведущий элемент = " << scientific << smallestPivot << endl;
        cerr << "Определитель матрицы = " << scientific << determinant << endl;
        cerr << "Решение может быть недостоверным." << endl;
    }

    vector<int> ipiv(dimension);
    int nrhs = 1;
    int lda = dimension;
    int ldb = dimension;
    int info = 0;

    auto startLapack = high_resolution_clock::now();
    dgesv_(&dimension, &nrhs, lapackMat.data(), &lda, ipiv.data(), lapackRhs.data(), &ldb, &info);
    auto endLapack = high_resolution_clock::now();
    auto lapackTime = duration_cast<microseconds>(endLapack - startLapack).count();

    if (info != 0) {
        cerr << "LAPACK: ошибка решения (info = " << info << ")" << endl;
        return 1;
    }

    vector<double> gaussRes = computeResidual(matrix, gaussSol, rhs);
    vector<double> lapackRes = computeResidual(matrix, lapackRhs, rhs);
    double normRhs = euclideanNorm(rhs);
    double gaussRelRes = euclideanNorm(gaussRes) / normRhs;
    double lapackRelRes = euclideanNorm(lapackRes) / normRhs;

    double errorNorm = 0.0;
    for (int i = 0; i < dimension; ++i) {
        double diff = gaussSol[i] - lapackRhs[i];
        errorNorm += diff * diff;
    }
    errorNorm = sqrt(errorNorm);

    cout << fixed << setprecision(6);
    cout << "Размерность системы: " << dimension << "\n";
    cout << "Время метода Гаусса: " << gaussTime << " мкс\n";
    cout << "Время LAPACK (dgesv): " << lapackTime << " мкс\n";
    cout << "Относительная невязка Гаусса: " << scientific << gaussRelRes << "\n";
    cout << "Относительная невязка LAPACK: " << scientific << lapackRelRes << "\n";
    cout << "Погрешность решения Гаусса относительно LAPACK (евклидова норма): " << scientific << errorNorm << "\n";

    return 0;
}