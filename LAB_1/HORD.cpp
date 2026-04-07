#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;


double F(double x) {
    return 2*x*x-x*x*x*x-1-log(x);
}

double df(double x) {
    return 4*x-4*x*x*x-1/x;
}

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

double F_round(double x, double delta) {
    return Round(F(x), delta);
}

double HORDA(double Left, double Right, double Eps, int &N) {
    double FLeft = F(Left);
    double FRight = F(Right);
    double X, Y;
    
    if (FLeft * FRight > 0.0) { 
        cout << "Неверное задание интервала\n"; 
        exit(1); 
    }
    if (Eps <= 0.0) { 
        cout << "Неверное задание точности\n"; 
        exit(1); 
    }
    
    N = 0;
    if (FLeft == 0.0) return Left;
    if (FRight == 0.0) return Right;
    
    do {
        X = Left - (Right - Left) * FLeft / (FRight - FLeft);
        Y = F(X);
        if (Y == 0.0) return X;
        
        if (Y * FLeft < 0.0) {
            Right = X;
            FRight = Y;
        } else {
            Left = X;
            FLeft = Y;
        }
        N++;
    } while ( fabs(Y) >= Eps);
    
    return X;
}


int main() {
    double a = 0.99, b = 1.15;
    double eps = 0.000001;

    // cout << "Введите левую границу a: ";
    // cin >> a;

    // cout << "Введите правую границу b: ";
    // cin >> b;

    // cout << "Введите точность ответа: ";
    // cin >> eps;


    int iterations;
    double root;
    
    cout << "МЕТОД ХОРД\n\n";
    cout << "Ответ = " << HORDA(a, b, eps, iterations) << endl;
    cout << "Количество итераций = " << iterations << endl << endl;

    
    cout << "ЧАСТЬ 1: ИССЛЕДОВАНИЕ СКОРОСТИ СХОДИМОСТИ\n";
    cout << "Eps\t\tКорень\t\tИтераций\t|f(x)|\n";
    cout << "------------------------------------------------\n";

    eps = 0.1;
    for (int i = 0; i < 10; i++) {
        root = HORDA(a, b, eps, iterations);
        cout << eps << "\t\t" << root << "\t\t" << iterations << "\t\t" << fabs(F(root)) << endl;
        eps /= 10;
    }
    




    cout << "\nЧАСТЬ 2: ИССЛЕДОВАНИЕ ОБУСЛОВЛЕННОСТИ\n";
    cout << "------------------------------------------------\n";
    double delta = eps;
    
    double exact_root = HORDA(a, b, 1e-10, iterations);
    double nu = 1.0 / fabs(df(exact_root));
    
    cout << "Точный корень: " << exact_root << endl;
    cout << "Число обусловленности ν_Δ = " << nu << endl;
    cout << "ν_Δ * Delta = " << nu * delta << endl;
    
    #define F(x) F_round(x, delta)
    
    double round_root = HORDA(a, b, 1e-9, iterations);
    double error = fabs(exact_root - round_root);
    
    cout << "Корень с округлением: " << round_root << endl;
    cout << "Разница |x_точн - x_округл| = " << error << endl;
    
    return 0;
}
