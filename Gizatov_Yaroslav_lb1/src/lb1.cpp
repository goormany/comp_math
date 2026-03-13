#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <windows.h>
#include <clocale>


using namespace std;


double f(double x) {
    return log(log(x)) - exp(-(x*x));
}


double Round(double val, double delta) {
    if (delta <= 0.0) return val;             
    return round(val / delta) * delta;          
}


double BISECT(double a, double b, double eps, int &iter) {
    double fa = f(a);
    double fb = f(b);

    
    if (fa * fb > 0) {
        cerr << "Ошибка: f(a) и f(b) одного знака – корень может отсутствовать." << endl;
        return 0.0;
    }

    iter = 0;
    double c;

    while ((b - a) >= 2.0 * eps) {
        c = (a + b) / 2.0;         
        double fc = f(c);

        if (fc == 0.0) break;

        if (fa * fc < 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
        ++iter;
    }

    return (a + b) / 2.0;
}


double BISECT_ROUND(double a, double b, double eps, double delta, int &iter) {
    double fa = Round(f(a), delta);
    double fb = Round(f(b), delta);

    if (fa * fb > 0) {
        cerr << "Ошибка: после округления f(a) и f(b) одного знака." << endl;
        return 0.0;
    }

    iter = 0;
    double c;

    while ((b - a) >= 2.0 * eps) {
        c = (a + b) / 2.0;
        double fc = Round(f(c), delta);

        if (fc == 0.0) break;

        if (fa * fc < 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
        ++iter;
    }

    return (a + b) / 2.0;
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    setlocale(LC_ALL, ".UTF-8");
    
    double a = 2.0, b = 3.0;

    double eps_base = 1e-10;


    int choice;
    cout << "1. Вычисление без округления\n";
    cout << "2. Зависимость числа итераций от точности Eps\n";
    cout << "3. Чувствительность к ошибкам округления\n";
    cout << "0. Выход\n";
    cout << "Ваш выбор: ";
    cin >> choice;   
    
    switch (choice) {
        
        case 1: {
            cout << "Введите точность eps (например, 1e-6): ";
            cin >> eps_base;

            int iter;
            double root = BISECT(a, b, eps_base, iter);

            cout << fixed << setprecision(10);
            cout << "Решение уравнения f(x)=0 методом бисекции" << endl;
            cout << "Интервал: [" << int(a) << ", " << int(b) << "]" << endl;
            cout << "Точность: " << eps_base << endl;
            cout << "Найденный корень: " << root << endl;
            cout << "Число итераций:   " << iter << "\n\n";
            break;
        }

        case 2: {
            cout << "Зависимость числа итераций от точности Eps" << endl;
            cout << setw(8) << "Eps" << setw(25) << "Итераций" << endl;
            cout << "---------------------------------------------------" << endl;


            for (double eps = 0.1; eps >= 1e-10; eps /= 10.0) {
                int it;
                double r = BISECT(a, b, eps, it);
                cout << setw(5) << eps << setw(10) << it << endl;
            }
            cout << endl;
            break;
        }

        case 3: {
            cout << "Чувствительность к ошибкам округления" << endl;
            cout << setw(11) << "Delta" << setw(32) << "Корень" << setw(25) << "Итераций" << endl;
            cout << "-------------------------------------------------------" << endl;

            vector<double> deltas = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 0.0};

            for (double delta : deltas) {
                int it;
                double r;

                if (delta == 0.0) {
                    r = BISECT(a, b, eps_base, it);
                } else {
                    r = BISECT_ROUND(a, b, eps_base, delta, it);
                }

                if (r == 0.0 && delta != 0.0) {
                    cout << setw(15) << delta << setw(25) << "ошибка" << setw(15) << "-" << endl;
                } else {
                    cout << setw(15) << delta << setw(25) << setprecision(10) << r
                        << setw(11) << it << endl;
                }
            }
            break;
        }
        default: break;
    }

    return 0;
}