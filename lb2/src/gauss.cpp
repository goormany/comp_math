#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iomanip>
using namespace std;


vector<double> gauss(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    
    // Прямой ход
    for (int i = 0; i < n; i++) {
        // Поиск главного элемента
        int maxRow = i;
        double maxVal = fabs(A[i][i]);
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > maxVal) {
                maxVal = fabs(A[k][i]);
                maxRow = k;
            }
        }
        
        if (maxVal < 1e-10) {
            throw runtime_error("Система вырождена");
        }
        
        // Перестановка строк
        if (maxRow != i) {
            swap(A[i], A[maxRow]);
            swap(b[i], b[maxRow]);
        }
        
        // Нормировка строки
        double div = A[i][i];
        for (int j = i; j < n; j++) {
            A[i][j] /= div;
        }
        b[i] /= div;
        
        // Вычитание из нижних строк
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    
    // Обратный ход
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
    }
    
    return x;
}

int main(int argc, char* argv[]) {
    string filename;
    string exp;
    // Получение имени файла из аргументов командной строки или запрос у пользователя
    if (argc > 1) {
        filename = argv[1];
        exp = argv[1];
    } else {
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
    
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);
    
    // Чтение матрицы
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (!(file >> A[i][j])) {
                cerr << "Ошибка: недостаточно данных в файле" << endl;
                return 1;
            }
        }
    }
    
    // Чтение вектора b
    for (int i = 0; i < n; i++) {
        if (!(file >> b[i])) {
            cerr << "Ошибка: недостаточно данных в файле" << endl;
            return 1;
        }
    }
    
    file.close();

    double eps = 1e-8;
    try {
        vector<double> x = gauss(A, b);

        double max_residual = 0;
        double sum_sq_residual = 0;
        double norm_b = 0;
        
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                sum += A[i][j] * x[j];
            }
            double residual = fabs(sum - b[i]);
            max_residual = max(max_residual, residual);
            sum_sq_residual += residual * residual;
            norm_b += b[i] * b[i];
        }
        
        double norm_l2 = sqrt(sum_sq_residual);
        double rel_residual = norm_l2 / sqrt(norm_b);
        
        cout << "\nСтатистика:" << endl;
        cout << "Невязка L2: " << scientific << norm_l2 << endl;
        cout << "Невязка L∞: " << max_residual << endl;
        cout << "Относительная невязка: " << rel_residual << endl;

        bool flag = true;
        
        // Сначала проверяем уравнения
        vector<double> errors(n);
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                sum += A[i][j] * x[j];
            }
            errors[i] = fabs(sum - b[i]);
            if (errors[i] > eps) {
                flag = false;
            }
        }
        
        cout << "\nРезультаты (первые 10):" << endl;
        for (int i = 0; i < min(n, 10); i++) {
            cout << "x" << i + 1 << " = " << setprecision(10) << x[i] 
                 << " ОШИБКА: " << errors[i] << endl;
        }
        
        cout << "\nПроверка:" << endl;
        if(flag){
            cout << "Проверка пройдена успешно!" << endl;
        }else{
            cout << "Проверка не пройдена!" << endl;
        }

        string out_file = (argc > 2) ? argv[2] : "gauss_solution.txt";
        ofstream fout(out_file);
        if (fout.is_open()) {
            fout << scientific << setprecision(12);
            fout << " L2 = " << norm_l2 << " Невязка = " << max_residual << endl;
            for (int i = 0; i < n; i++)
                fout << "x" << i+1 << " = " << x[i] << endl;
            fout.close();
            cout << "\nСохранено в " << out_file << endl;
        }
            
    } catch (const exception& e) {
        cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}

// g++ gauss.cpp -o gauss.exe && ./gauss.exe 