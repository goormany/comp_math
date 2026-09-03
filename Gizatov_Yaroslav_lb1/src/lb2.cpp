#include <iostream>
#include <cmath>
#include <iomanip>
#include <windows.h>
#include <clocale>

using namespace std;


double Round(double x, double delta) {
    if (delta <= 0.0) return x;        
    return round(x / delta) * delta;
}


double f(double x) {
    return log(log(x)) - exp(-(x*x));
}

double f_round(double x, double delta) {
    return Round(f(x), delta);
}

bool HORDA(double a, double b, double eps, double delta, int &iter, double &c) {
    const int MAX_ITER = 1000000;               
    iter = 0;

    double fa, fb;
    
    if (delta > 0) {
        fa = f_round(a, delta);
        fb = f_round(b, delta);
    } else {
        fa = f(a);
        fb = f(b);
    }

    if (fa * fb >= 0.0) {
        cout << "Ошибка: f(a) и f(b) должны иметь противоположные знаки." << endl;
        cout << "f(a) = " << fa << ", f(b) = " << fb << endl;
        return false;
    }

    do {
        c = a - fa * (b - a) / (fb - fa);

        double fc;
        if (delta > 0) {
            fc = f_round(c, delta);
        } else {
            fc = f(c);
        }
        iter++;

        if (fabs(fc) < eps) {
            return true;
        }

        if (fc > 0.0) {
            a = c;
            fa = fc;
        } else {
            b = c;
            fb = fc;
        }

        if (iter >= MAX_ITER) {
            cout << "Достигнуто максимальное число итераций." << endl;
            return false;
        }
    } while (true);
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    setlocale(LC_ALL, ".UTF-8");

    double a = 2;            
    double b = 3;           
    double eps = 1e-6;          
    double delta = 1e-8;

    int iter;
    double root;


    int choice;
    cout << "1. Вычисление без округления\n";
    cout << "2. Вычисление с округлением\n";
    cout << "3. Исследование скорости сходимости (без округления)\n";
    cout << "4. Исследование обусловленности (влияние точности округления)\n";
    cout << "0. Выход\n";
    cout << "Ваш выбор: ";
    cin >> choice;   
    
    switch (choice) {

        case 1: {
            
            cout << "Введите точность eps (например, 1e-6): ";
            cin >> eps;    

            cout << fixed << setprecision(10);
            cout << "Начальный отрезок: [" << int(a) << ", " << int(b) << "]\n";
            cout << "Точность eps = " << eps << "\n\n";
            
            cout << "Вычисление без округления (delta = 0)\n";
            if (HORDA(a, b, eps, 0.0, iter, root)) {
                cout << "Корень найден: x = " << root << endl;
                cout << "Значение функции f(x) = " << f(root) << endl;
                cout << "Число итераций: " << iter << "\n\n";
            }
            break;
        }
        case 2: {

            cout << "Введите точность eps (например, 1e-6): ";
            cin >> eps;

            cout << "Вычисление с округлением (delta = " << delta << ")\n";
            if (HORDA(a, b, eps, delta, iter, root)) {
                cout << "Корень найден: x = " << root << endl;
                cout << "Значение функции f_raw(x) = " << f(root) << endl;
                cout << "Значение функции f_round(x) = " << f_round(root, delta) << endl;
                cout << "Число итераций: " << iter << "\n\n";
            }
            break;
        }

        case 3: {
            cout << "Исследование скорости сходимости (без округления)\n";
            double eps_vals[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
            for (double e : eps_vals) {
                int it;
                double r;
                if (HORDA(a, b, e, 0.0, it, r)) {
                    cout << "eps = " << e << ": итераций = " << it
                        << ", корень = " << r << endl;
                }
            }
            cout << endl;
            break;
        }

        case 4: {
            cout << "Исследование обусловленности (влияние точности округления)\n";
            double delta_vals[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
            for (double d : delta_vals) {
                int it;
                double r;
                if (HORDA(a, b, eps, d, it, r)) {
                    cout << "delta = " << d << ": итераций = " << it
                        << ", корень = " << r
                        << ", f_raw(x) = " << f(r)
                        << ", f_round(x) = " << f_round(r, d) << endl;
                } else {
                    cout << "delta = " << d << ": не удалось найти корень" << endl;
                }
            }
            break;
        }
        default: break;
    }
    return 0;
}