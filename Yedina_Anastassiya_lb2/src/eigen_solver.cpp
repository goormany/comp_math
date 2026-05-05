#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Использование: " << argv[0] << " <файл_с_матрицей>\n";
        return 1;
    }

    ifstream in(argv[1]);
    if (!in.is_open()) {
        cerr << "ошибка: не удалось открыть файл.\n";
        return 1;
    }

    int n;
    in >> n;

    MatrixXd A(n, n);
    VectorXd b(n);


    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            in >> A(i, j);
        }
        in >> b(i);
    }
    in.close();

    cout << "размер матрицы: " << n << "x" << n << "\n";
    cout << "решаем систему методом Eigen\n";

    auto start_time = chrono::high_resolution_clock::now();

    VectorXd x = A.partialPivLu().solve(b);

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;

    cout << "Eigen решил систему за " << elapsed.count() << " сек.\n";

    VectorXd residual = A * x - b;
    double norm_L2 = residual.norm();
    double norm_Linf = residual.lpNorm<Infinity>();

    cout << "\nрезультаты EIGEN\n";
    cout << scientific << setprecision(6);
    cout << "норма невязки L2: " << norm_L2 << "\n";
    cout << "норма невязки L_inf: " << norm_Linf << "\n";

    return 0;
}