#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
using namespace std;


double F(double x);
double F1(double x);
double F2(double x);
double NEWTON (double X,double Eps,int &N);
double Round(double X, double Delta);
double delta = 1E-8;

double F(double x){
    return Round(tan(x) - 1/x, delta);
}

double F1(double x){
    return Round(1/(cos(x)*cos(x)) + 1/(x*x), delta);
}

double F2(double x){
    double a = x*x*x;
    double b = cos(x)*cos(x)*cos(x);
    return Round(2*(a*sin(x) - b)/(a*b), delta);
}

double NEWTON (double X,double Eps,int &N)
{
    extern double F1 (double);
    double Y,Y1,DX;
    N=0;
    do{
        Y  = F(X);
        if (Y==0.0) return (X);
        Y1 = F1(X);
        if (Y1==0.0) {puts("Производная обратилась в ноль\n");exit(1);}
        DX=Y/Y1; X=X-DX; N++;
    }
    while (fabs(DX)>Eps);
    return (X);
}


double Round(double X, double Delta){
    if (Delta <= 1E-9) { 
        puts("Неверное задание точности округления\n"); 
        exit(1); 
    }
    if (X > 0.0) 
        return (Delta * (long((X / Delta) + 0.5)));
    else    
        return (Delta * (long((X / Delta) - 0.5)));
}

int main(){
    double Left = 0.5;
    double Right = 1.0;
    double Eps;
    double Root;
    int N;

    cout << "Уравнение: tg(x) - 1/x = 0" << std::endl;
    cout << "Интервал поиска: [" << Left << ", " << Right << "]" << std::endl;
    cout << std::endl;

    double x;
    cout << "Введите начальное приближение: ";
    cin >> x;
    
    std::cout << "Производная в точке " << x << ": " << F1(x) << std::endl;
    std::cout << "Вторая производная в точке " << x << ": " << F2(x) << std::endl;

    if(F(x)*F2(x) > 0){
        cout << "Хорошее начальное приближение, выполнено условие сходимости" << endl;
    }else{
        cout << "Плохое начальное приближение, условие сходимости не выполнено" << endl;
        return 1;
    }
    cout << F(x)*F2(x) << endl;
    printf("--------------------------------------------------------\n");
    printf("    Eps      | Итерации |      Корень      | Невязка\n");
    printf("--------------------------------------------------------\n");
    
    Eps = 0.1;
    while (Eps >= 0.00000001){
        double Root = NEWTON(x, Eps, N);
        double Residual = fabs(F(Root));
        printf("%10.7f | %8d | %16.12f | %10.3e\n", Eps, N, Root, Residual);
        Eps /= 10.0;
    }

    printf("\n\nИсследование чувствительности к ошибкам в исходных данных:\n");
    printf("Точность метода Eps = 1e-8\n");
    printf("----------------------------------------------------------------\n");
    printf("   Delta    | Итерации |      Корень      | Ошибка\n");
    printf("----------------------------------------------------------------\n");
    double TrueEps = 1e-8;
    int TrueNit;
    double TrueRoot = NEWTON(x, TrueEps, TrueNit);

    for (double d = 0.1; d >= 1e-8; d /= 10) {
        delta = d;
        int Nit = 0;
        double root = NEWTON(x, TrueEps, Nit);  
        double Error = fabs(root - TrueRoot);
        printf("%10.7f | %8d | %16.12f | %10.3e\n", 
            delta, Nit, root, Error);
    }
    
    return 0;
}