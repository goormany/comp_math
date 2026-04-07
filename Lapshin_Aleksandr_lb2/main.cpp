#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm> 
#include <Eigen/Dense>
#include <chrono>

using namespace std;
const float EPS = 1e-9f;
int n;



vector<float> method_gaus(vector<vector<float>> A, vector<float> b) {
    vector<float> x(n);
    for (int i = 0; i < n; i++) {
        float current_max = 0.0f;
        int pivot_row = -1;

        for (int k = i; k < n; k++) {
            float val = fabs(A[k][i]);
            if (val > current_max) {
                current_max = val;
                pivot_row = k;
            }
        }

        if (current_max < EPS) {
            cout << "Матрица вырождена" << endl;
            return vector<float>();
        }
        
        if (pivot_row != i) {
            swap(A[i], A[pivot_row]);
            swap(b[i], b[pivot_row]);
        }

        for (int k = i + 1; k < n; k++) {
            float ratio = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= A[i][j] * ratio;
            }
            b[k] -= b[i] * ratio;
        }
    }
    
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--) {
        float sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= A[i][j] * x[j];
        }
        x[i] = sum / A[i][i];
    }
    
    return x;
}



void print_result(vector<float> result_gauss, chrono::nanoseconds time_result_gaus, Eigen::VectorXf result_eigen, chrono::nanoseconds time_result_eigen) {
    cout << "\nРЕЗУЛЬТАТЫ" << endl;
    cout << "Метод Гаусса (частичный выбор):" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << result_gauss[i] << endl;
    }
    
    cout << "\nМетод встроенной библиотеки (Eigen):" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << result_eigen(i) << endl;
    }
    
    float L1_norm = 0.0f;
    float L2_norm = 0.0f;
    float Linf_norm = 0.0f;
    
    for (int i = 0; i < n; i++) {
        float diff = fabs(result_gauss[i] - result_eigen(i));
        L1_norm += diff;
        L2_norm += diff * diff;
        if (diff > Linf_norm) {
            Linf_norm = diff;
        }
    }

    L2_norm = sqrt(L2_norm);
    
    cout << "\nСравнение норм:" << endl;
    cout << "L1-норма (Манхэттенское расстояние: L1 = |x₁ - x₂| + |y₁ - y₂| + |z₁ - z₂| + |w₁ - w₂|) = " << L1_norm << endl;
    cout << "L2-норма (Евклидово расстояние: L2 = √[(x₁-x₂)² + (y₁-y₂)² + (z₁-z₂)² + (w₁-w₂)²]) = " << L2_norm << endl;
    cout << "L∞-норма (Самое большое отклонение: L∞ = max(|x₁-x₂|, |y₁-y₂|, |z₁-z₂|, |w₁-w₂|)) = " << Linf_norm << endl;
    
    cout << "\nСравнение времени:" << endl;
    cout << "Метод Гаусса: " << time_result_gaus.count() << " ns" << endl;
    cout << "Метод Eigen: " << time_result_eigen.count() << " ns" << endl;
}



int main() {
    cout << "Решение уравненй вида Ax=b\n";
    cout << "Введите размер матрицы A: ";
    cin >> n;
    
    cout << "Введите матрицу A (" << n << " x " << n << "):" << endl;
    vector<vector<float>> A(n, vector<float>(n));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    vector<float> b(n);
    cout << "Введите вектор правой части b (" << n << " элементов):" << endl;
    for (int i = 0; i < n; i++) {
        cin >> b[i];
    }
    

    auto start_gauss = chrono::high_resolution_clock::now();
    vector<float> result_gauss = method_gaus(A, b);
    auto end_gauss = chrono::high_resolution_clock::now();
    chrono::nanoseconds time_result_gauss = chrono::duration_cast<chrono::nanoseconds>(end_gauss - start_gauss);
    

    if (result_gauss.empty()) {
        cout << "Ошибка в методе Гаусса!" << endl;
        return 1;
    }
    
    Eigen::MatrixXf A_eigen(n, n);
    Eigen::VectorXf b_eigen(n);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A_eigen(i, j) = A[i][j];
        }
        b_eigen(i) = b[i];
    }
    


    auto start_eigen = chrono::high_resolution_clock::now();
    Eigen::VectorXf result_eigen = A_eigen.colPivHouseholderQr().solve(b_eigen);
    auto end_eigen = chrono::high_resolution_clock::now();
    chrono::nanoseconds time_result_eigen = chrono::duration_cast<chrono::nanoseconds>(end_eigen - start_eigen);
    


    print_result(result_gauss, time_result_gauss, result_eigen, time_result_eigen);
    return 0;
}

