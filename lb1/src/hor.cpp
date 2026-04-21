#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>
using namespace std;


double delta = 0.0;

double Round(double X, double Delta) {
    if (Delta <= 1E-9) {
        cout << "Неверное задание точности округления" << std::endl;
        exit(1);
    }
    if (X > 0.0)
        return (Delta * (long((X / Delta) + 0.5)));
    else
        return (Delta * (long((X / Delta) - 0.5)));
}

double F(double x) {
    double fx = tan(x) - 1.0 / x;
    if (delta == 0.0) {
        return fx;
    }
    return Round(fx, delta);
}

double HORDA(double Left, double Right, double Eps, int &N) {
    double FLeft = F(Left);
    double FRight = F(Right);
    double X, Y;

    if (FLeft * FRight > 0.0) {
        cout << "Неверное задание интервала" << std::endl;
        exit(1);
    }
    if (Eps <= 0.0) {
        cout << "Неверное задание точности" << std::endl;
        exit(1);
    }

    N = 0;
    if (FLeft == 0.0) return Left;
    if (FRight == 0.0) return Right;

    do {
        X = Left - (Right - Left) * FLeft / (FRight - FLeft);
        Y = F(X);
        N++;
        if (Y == 0.0) return X;
        if (Y * FLeft < 0.0) {
            Right = X;
            FRight = Y;
        } else {
            Left = X;
            FLeft = Y;
        }
    } while (fabs(Y) >= Eps);

    return X;
}

int main() {
    double left = 0.5;
    double right = 1.0;
    int N;

    cout << "Уравнение: tg(x) - 1/x = 0" << std::endl;
    cout << "Интервал поиска: [" << left << ", " << right << "]" << std::endl;
    cout << std::endl;

    if (F(left) * F(right) > 0) {
        cout << "На интервале [" << left << ", " << right
                  << "] нет гарантии наличия корня!" << std::endl;
        cout << "F(left) = " << F(left) << ", F(right) = " << F(right) << std::endl;
        return 1;
    }

    printf("Пункт 4. Зависимость числа итераций от точности Eps:\n");
    printf("--------------------------------------------------------\n");
    printf("    Eps      | Итерации |      Корень      | Невязка\n");
    printf("--------------------------------------------------------\n");
    
    double Eps = 0.1;
    while (Eps >= 0.00000001){
        double Root = HORDA(left, right, Eps, N);
        double Residual = fabs(F(Root));
        printf("%10.7f | %8d | %16.12f | %10.3e\n", Eps, N, Root, Residual);
        Eps /= 10.0;
    }

    printf("\n\nПункт 5. Исследование чувствительности к ошибкам в исходных данных:\n");
    printf("Точность метода Eps = 1e-8\n");
    printf("----------------------------------------------------------------\n");
    printf("   Delta    | Итерации |      Корень      | Ошибка\n");
    printf("----------------------------------------------------------------\n");
    double TrueEps = 1e-8;
    int TrueNit;
    double TrueRoot = HORDA(left, right, TrueEps, TrueNit);

    for (double d = 0.1; d >= 0.000001; d /= 10) {
        delta = d;
        int Nit = 0;
        double root = HORDA(left, right, TrueEps, Nit);  
        double Error = fabs(root - TrueRoot);
        printf("%10.7f | %8d | %16.12f | %10.3e\n", 
            delta, Nit, root, Error);
    }
    return 0;
}