#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <chrono>
#include <iomanip>

using namespace std;

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

    Eigen::MatrixXd A(n, n);
    Eigen::VectorXd b(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(fin >> A(i, j))) {
                cerr << "Ошибка: недостаточно данных." << endl;
                return 1;
            }
        }
        if (!(fin >> b(i))) {
            cerr << "Ошибка: недостаточно данных." << endl;
            return 1;
        }
    }
    fin.close();

    auto start = chrono::high_resolution_clock::now();
    Eigen::VectorXd x = A.lu().solve(b);
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double>(end - start).count();

    Eigen::VectorXd residual = A * x - b;
    double norm_l2 = residual.norm();          // L2 норма
    double norm_linf = residual.lpNorm<Eigen::Infinity>(); // L∞ норма

    ofstream fout(output_file);
    if (!fout.is_open()) {
        cerr << "Ошибка: не удалось создать выходной файл " << output_file << endl;
        return 1;
    }

    fout << fixed << setprecision(12);
    for (int i = 0; i < n; ++i) {
        fout << "x" << i + 1 << " = " << x(i) << endl;
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