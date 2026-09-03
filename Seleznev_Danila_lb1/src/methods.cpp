#include <stdio.h>
#include <math.h>
#include <stdlib.h>

extern double F(double);
extern double F1(double);
extern double F2(double);
extern double Phi(double);

// Округление значения X с точностью Delta
double Round(double X, double Delta)
{
    if (Delta <= 1E-9) {
        puts("Неверное задание точности округления\n");
        exit(1);
    }
    if (X > 0.0)
        return (Delta * (long((X / Delta) + 0.5)));
    else
        return (Delta * (long((X / Delta) - 0.5)));
}

// Метод деления отрезка пополам (бисекция)
double BISECT(double Left, double Right, double Eps, int &N)
{
    double E = fabs(Eps) * 2.0;
    double FLeft = F(Left);
    double FRight = F(Right);
    double X = (Left + Right) / 2.0;
    double Y;

    if (FLeft * FRight > 0.0) {
        puts("Неверное задание интервала\n");
        exit(1);
    }
    if (Eps <= 0.0) {
        puts("Неверное задание точности\n");
        exit(1);
    }

    N = 0;
    if (FLeft == 0.0)  return Left;
    if (FRight == 0.0) return Right;

    while ((Right - Left) >= E)
    {
        X = 0.5 * (Right + Left);
        Y = F(X);
        if (Y == 0.0) return X;
        
        if (Y * FLeft < 0.0)
            Right = X;
        else {
            Left = X;
            FLeft = Y;
        }
        N++;
    }
    return X;
}

// Метод хорд
double HORDA(double Left, double Right, double Eps, int &N)
{
    double FLeft = F(Left);
    double FRight = F(Right);
    double X, Y;

    if (FLeft * FRight > 0.0) {
        puts("Неверное задание интервала\n");
        exit(1);
    }
    if (Eps <= 0.0) {
        puts("Неверное задание точности\n");
        exit(1);
    }

    N = 0;
    if (FLeft == 0.0)  return Left;
    if (FRight == 0.0) return Right;

    do
    {
        X = Left - (Right - Left) * FLeft / (FRight - FLeft);
        Y = F(X);
        if (Y == 0.0) return X;
        
        if (Y * FLeft < 0.0) {
            Right = X;
            FRight = Y;
        }
        else {
            Left = X;
            FLeft = Y;
        }
        N++;
    }
    while (fabs(Y) >= Eps);

    return X;
}

// Метод Ньютона (касательных)
double NEWTON(double X, double Eps, int &N, double Left, double Right)
{
    double Y, Y1, DX;
    double XPrev;
    N = 0;
    
    // Вычисление m1 и M2 для улучшенного критерия остановки
    double m1 = 1e10, M2 = 0.0;
    for (double t = Left; t <= Right; t += 0.001) {
        double f1 = fabs(F1(t));
        double f2 = fabs(F2(t));
        if (f1 < m1) m1 = f1;
        if (f2 > M2) M2 = f2;
    }
    
    double Eps0 = sqrt(2.0 * m1 * Eps / M2);
    
    do
    {
        Y = F(X);
        if (Y == 0.0) return X;

        Y1 = F1(X);
        if (Y1 == 0.0) {
            puts("Производная обратилась в ноль\n");
            exit(1);
        }

        XPrev = X;
        DX = Y / Y1;
        X = X - DX;
        N++;
        
        // Улучшенный критерий остановки
        if (N > 1 && fabs(X - XPrev) < Eps0) {
            return X;
        }
    }
    while (fabs(DX) > Eps);
    
    return X;
}

// Метод простых итераций
double ITER(double X0, double Eps, int &N)
{
    if (Eps <= 0.0) {
        puts("Неверное задание точности\n");
        exit(1);
    }
    
    double XPrev, Y;
    N = 0;
    
    do
    {
        XPrev = X0;
        Y = Phi(XPrev);
        X0 = Y;
        N++;
    }
    while (fabs(Y - XPrev) >= Eps);
    
    return Y;
}