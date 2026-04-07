#include <iostream>
#include <cmath>
using namespace std;


double Round(double X, double Delta) {
    if (Delta <= 1E-9) { 
        cout << "Неверное задание точности округления\n"; 
        exit(1); 
    }
    if (X > 0.0) 
        return (Delta * (long((X/Delta) + 0.5)));
    else 
        return (Delta * (long((X/Delta) - 0.5)));
}

double F(double x) {
    return 2*x*x - x*x*x*x - 1 - log(x);
}

double dF(double x) {
    return 4*x - 4*x*x*x - 1/x;
}

double ddF(double x) {
    return 4 - 12*x*x + 1/(x * x);
}

double F_round(double x, double delta) {
    return Round(F(x), delta);
}

double NEWTON(double X, double Eps, int &N) {
    double DX = F(X) / dF(X);
    N = 0;
    do {
        double Y = F(X);
        if (Y == 0.0) return X;
        double Y1 = dF(X);
        if (Y1 == 0.0) {
            cout << "Производная обратилась в ноль\n";
            exit(1);
        }
        DX = Y / Y1;
        X = X - DX;
        N++;
    } while (fabs(DX) > Eps);
    return X;
}

int main() {
    double x0 = 1.15;
    // double x0 = 0.55;

    cout << F(x0) * ddF(x0) << endl;

    double eps = 1e-6;
    double Delta = 0.000001;
    int iterations;
    double root;

    // cout << "Введите начальное значение x0: ";
    // cin >> x0;

    // cout << "Введите точность ответа: ";
    // cin >> eps;

    cout << "МЕТОД НЬЮТОНА\n";
    cout << "Ответ = " << NEWTON(x0, eps, iterations) << endl;
    cout << "Количество итераций = " << iterations << endl << endl;


    cout << "ЧАСТЬ 1: СКОРОСТЬ СХОДИМОСТИ\n";
    cout << "Eps\t\tКорень\t\tИтераций\t|F(x)|\n";
    cout << "------------------------------------------------\n";

    eps = 0.1;
    for (int i = 0; i < 10; i++) {
        root = NEWTON(x0, eps, iterations);
        cout << eps << "\t\t" << root << "\t\t" << iterations << "\t\t" << fabs(F(root)) << endl;
        eps /= 10;
    }



    cout << "\nЧАСТЬ 2: ИССЛЕДОВАНИЕ ОБУСЛОВЛЕННОСТИ\n";
    cout << "------------------------------------------------\n";
    double delta = eps;

    double exact_root = NEWTON(x0, 1e-10, iterations);
    double nu = 1.0 / fabs(dF(exact_root));
    
    cout << "Точный корень: " << exact_root << endl;
    cout << "Число обусловленности ν_Δ = " << nu << endl;
    cout << "ν_Δ * Delta = " << nu * delta << endl;

    #define F(x) F_round(x, delta)
    
    double round_root = NEWTON(x0, 1e-9, iterations);
    double error = fabs(exact_root - round_root);
    
    cout << "Корень с округлением: " << round_root << endl;
    cout << "Разница |x_точн - x_округл| = " << error << endl;
    
    return 0;
}
