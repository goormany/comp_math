#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <random>

extern "C" {
    void dgesv_(int* n, int* nrhs, double* A, int* lda, int* ipiv,
                double* B, int* ldb, int* info);
}

using namespace std;
using namespace chrono;

vector<vector<double>> generateDiagDominantMatrix(int dim) {
    unsigned seed = system_clock::now().time_since_epoch().count();
    mt19937 rng(seed);
    uniform_real_distribution<double> uniformDist(-1.0, 1.0);
    uniform_real_distribution<double> positiveDist(0.1, 1.0);

    vector<vector<double>> mat(dim, vector<double>(dim));
    for (int i = 0; i < dim; ++i) {
        double offDiagSum = 0.0;
        for (int j = 0; j < dim; ++j) {
            if (j != i) {
                mat[i][j] = uniformDist(rng);
                offDiagSum += fabs(mat[i][j]);
            }
        }
        mat[i][i] = offDiagSum + positiveDist(rng);
    }
    return mat;
}

vector<double> computeRHS(const vector<vector<double>>& mat) {
    int dim = mat.size();
    vector<double> rhs(dim, 0.0);
    for (int i = 0; i < dim; ++i) {
        double sum = 0.0;
        for (int j = 0; j < dim; ++j)
            sum += mat[i][j] * (j + 1);
        rhs[i] = sum;
    }
    return rhs;
}

double vectorNorm(const vector<double>& vec) {
    double sqSum = 0.0;
    for (double val : vec) sqSum += val * val;
    return sqrt(sqSum);
}

vector<double> computeResidual(const vector<vector<double>>& mat,
                               const vector<double>& sol,
                               const vector<double>& rhs) {
    int dim = mat.size();
    vector<double> residual(dim, 0.0);
    for (int i = 0; i < dim; ++i) {
        double sum = 0.0;
        for (int j = 0; j < dim; ++j)
            sum += mat[i][j] * sol[j];
        residual[i] = rhs[i] - sum;
    }
    return residual;
}

bool isDiagonallyDominant(const vector<vector<double>>& mat) {
    int dim = mat.size();
    for (int i = 0; i < dim; ++i) {
        double diag = fabs(mat[i][i]);
        double sum = 0.0;
        for (int j = 0; j < dim; ++j)
            if (j != i) sum += fabs(mat[i][j]);
        if (diag <= sum) return false;
    }
    return true;
}

vector<double> solveJacobi(const vector<vector<double>>& mat,
                           const vector<double>& rhs,
                           double tolerance, int maxIters,
                           int& iterationsDone, bool& success) {
    int dim = mat.size();
    vector<double> prevX(dim, 0.0);
    vector<double> nextX(dim);

    for (int i = 0; i < dim; ++i) {
        if (fabs(mat[i][i]) < 1e-12) {
            success = false;
            return {};
        }
    }

    for (iterationsDone = 1; iterationsDone <= maxIters; ++iterationsDone) {
        for (int i = 0; i < dim; ++i) {
            double sum = 0.0;
            for (int j = 0; j < dim; ++j)
                if (j != i) sum += mat[i][j] * prevX[j];
            nextX[i] = (rhs[i] - sum) / mat[i][i];
        }

        double maxChange = 0.0;
        for (int i = 0; i < dim; ++i) {
            double change = fabs(nextX[i] - prevX[i]);
            if (change > maxChange) maxChange = change;
        }
        if (maxChange < tolerance) {
            success = true;
            return nextX;
        }
        prevX = nextX;
    }
    success = false;
    return nextX;
}

int main() {
    int dim;
    cout << "Введите размер матрицы n: ";
    cin >> dim;
    if (dim <= 0) {
        cerr << "Ошибка: размер должен быть положительным." << endl;
        return 1;
    }

    vector<vector<double>> matrix = generateDiagDominantMatrix(dim);
    vector<double> rhs = computeRHS(matrix);

    vector<double> exactSolution(dim);
    for (int j = 0; j < dim; ++j)
        exactSolution[j] = j + 1;

    cout << "Сгенерирована система " << dim << "x" << dim;
    if (isDiagonallyDominant(matrix))
        cout << " (строгое диагональное преобладание).\n";
    else
        cout << " (диагональное преобладание отсутствует).\n";

    const double TOLERANCE = 1e-3;
    const int MAX_ITERATIONS = 50000;
    int iterations = 0;
    bool convergedFlag = false;

    auto startJacobi = high_resolution_clock::now();
    vector<double> jacobiSolution = solveJacobi(matrix, rhs, TOLERANCE, MAX_ITERATIONS,
                                                 iterations, convergedFlag);
    auto endJacobi = high_resolution_clock::now();
    auto jacobiTime = duration_cast<microseconds>(endJacobi - startJacobi).count();

    vector<double> lapackMatrix(dim * dim);
    for (int j = 0; j < dim; ++j)
        for (int i = 0; i < dim; ++i)
            lapackMatrix[j * dim + i] = matrix[i][j];

    vector<double> lapackSolution = rhs;
    vector<int> pivotIndices(dim);
    int numRHS = 1;
    int lda = dim, ldb = dim;
    int info = 0;

    auto startLapack = high_resolution_clock::now();
    dgesv_(&dim, &numRHS, lapackMatrix.data(), &lda, pivotIndices.data(),
           lapackSolution.data(), &ldb, &info);
    auto endLapack = high_resolution_clock::now();
    auto lapackTime = duration_cast<microseconds>(endLapack - startLapack).count();

    if (info != 0) {
        cerr << "Ошибка LAPACK dgesv: info = " << info << "\n";
        return 1;
    }

    cout << fixed << setprecision(6);
    cout << "\nРазмерность системы: " << dim << "\n";
    cout << "Время метода Якоби: " << jacobiTime << " мкс";
    if (!convergedFlag) cout << " (НЕ сошёлся)";
    cout << "\n";
    cout << "Время LAPACK (dgesv): " << lapackTime << " мкс\n";

    if (convergedFlag)
        cout << "Метод Якоби сошёлся за " << iterations << " итераций\n";
    else
        cout << "Метод Якоби НЕ сошёлся за " << MAX_ITERATIONS << " итераций\n";

    double normRHS = vectorNorm(rhs);

    if (!jacobiSolution.empty()) {
        vector<double> residualJacobi = computeResidual(matrix, jacobiSolution, rhs);
        cout << "Относительная невязка Якоби: " << scientific
             << vectorNorm(residualJacobi) / normRHS << "\n";
    }

    vector<double> residualLapack = computeResidual(matrix, lapackSolution, rhs);
    cout << "Относительная невязка LAPACK: " << scientific
         << vectorNorm(residualLapack) / normRHS << "\n";

    cout << "\nТочное решение: x_j = j+1 (j = 0.." << dim-1 << ")\n";

    if (!jacobiSolution.empty()) {
        double errorJacobi = 0.0;
        for (int i = 0; i < dim; ++i)
            errorJacobi += (jacobiSolution[i] - exactSolution[i]) *
                           (jacobiSolution[i] - exactSolution[i]);
        errorJacobi = sqrt(errorJacobi);
        cout << "Погрешность ||x_jacobi - x_exact||_2 : " << scientific << errorJacobi << "\n";
    }

    double errorLapack = 0.0;
    for (int i = 0; i < dim; ++i)
        errorLapack += (lapackSolution[i] - exactSolution[i]) *
                       (lapackSolution[i] - exactSolution[i]);
    errorLapack = sqrt(errorLapack);
    cout << "Погрешность ||x_lapack - x_exact||_2 : " << scientific << errorLapack << "\n";

    return 0;
}