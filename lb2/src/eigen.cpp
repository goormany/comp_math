#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <Eigen/SVD>
#include <Eigen/Dense>

using namespace std;
using namespace chrono;

int main(int argc, char* argv[]) {
    string filename = (argc > 1) ? argv[1] : "";
    if (filename.empty()) {
        cout << "Введите имя файла: ";
        cin >> filename;
    }
    
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Ошибка: не удалось открыть файл " << filename << endl;
        return 1;
    }
    
    int n;
    file >> n;
    if (n <= 0) {
        cerr << "Ошибка: неверный размер матрицы" << endl;
        return 1;
    }
    
    Eigen::MatrixXd A(n, n);
    Eigen::VectorXd b(n);
    
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            file >> A(i, j);
    
    for (int i = 0; i < n; i++)
        file >> b(i);
    
    file.close();
    
    auto start = high_resolution_clock::now();
    
    // SVD решение
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
    
    // Адаптивный порог
    double threshold = (cond > 1e12) ? 1e-10 : 1e-12;
    svd.setThreshold(threshold);
    
    Eigen::VectorXd x = svd.solve(b);
    
    auto end = high_resolution_clock::now();
    double elapsed = duration<double>(end - start).count();
    
    // Невязки
    Eigen::VectorXd residual = A * x - b;
    double norm_l2 = residual.norm();
    double norm_linf = residual.lpNorm<Eigen::Infinity>();
    double rel_residual = norm_l2 / b.norm();
    
    // Вывод
    cout << fixed << setprecision(10);
    cout << "\nРешение (первые 10):" << endl;
    for (int i = 0; i < min(n, 10); i++)
        cout << "x" << i+1 << " = " << x(i) << endl;
    
    cout << "\nСтатистика:" << endl;
    cout << "Число обусловленности: " << scientific << cond << endl;
    cout << "Невязка L2: " << norm_l2 << endl;
    cout << "Относительная невязка: " << rel_residual << endl;
    cout << "Время: " << fixed << setprecision(3) << elapsed << " сек." << endl;
    
    string out_file = (argc > 2) ? argv[2] : "eigen_solution.txt";
    ofstream fout(out_file);
    if (fout.is_open()) {
        fout << scientific << setprecision(12);
        fout << "# cond = " << cond << " L2 = " << norm_l2 << " время = " << elapsed << endl;
        for (int i = 0; i < n; i++)
            fout << "x" << i+1 << " = " << x(i) << endl;
        fout.close();
        cout << "\nСохранено в " << out_file << endl;
    }
    
    return 0;
}

// g++ -std=c++17 -I /mnt/c/msys/ucrt64/include/eigen3 eigen.cpp -o eigen.exe