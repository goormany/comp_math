#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;


double F(double x){
    return tan(x) - 1/x;
}

double F1(double x){
    return 1/(cos(x)*cos(x)) + 1/(x*x);
}

double Round(double X, double Delta) {
    if (Delta <= 1e-12) { 
        puts("Неверное задание точности округления\n"); 
    }
    if (X > 0.0) 
        return (Delta * (long((X / Delta) + 0.5)));
    else    
        return (Delta * (long((X / Delta) - 0.5)));
}

double phi(double x) {
    return x-0.2056*(tan(x)-1/x);
}

double phi_deriv(double x) {
    return 1-0.2056*(1/(cos(x)*cos(x))+1/(x*x));
}

double ITER(double X0, double Eps, double delta, int &N){
  if (Eps<=0.0) {puts("Неверное задание точности\n");exit (1);}
  double X1=Round(phi(X0), delta);
  double X2=Round(phi(X1), delta);
  N = 2;
  int max_iter = 100;
  while( fabs(X1 - X2) > Eps && N < max_iter){
	X0 = X1;
	X1 = X2;
	X2 = Round(phi(X1), delta);
        N++;
  }
   return(X2);
}


int main() {
    int g = rand();
    cout << g << endl;
    double Delta = 1e-8;
    cout << "Уравнение: tg(x) - 1/x = 0\n";
    cout << "Интервал с корнем: [0.5, 1.0]\n\n";

    double x_test1 = 0.5;
    double x_test2 = 1;
    
    cout << F1(0.5) << "  " << F1(1) << endl;
    cout << "φ₁(x) = arctg(1/x)\n";
    cout << "φ₁(0.5) = " << phi(0.5) << endl;
    cout << "φ₁'(0.5) = " << phi_deriv(0.5) << endl;
    cout << "|φ₁'(0.5)| = " << fabs(phi_deriv(0.5)) << endl;

    cout << "φ₁(1) = " << phi(1) << endl;
    cout << "φ₁'(1) = " << phi_deriv(1) << endl;
    cout << "|φ₁'(1)| = " << fabs(phi_deriv(1)) << endl;
    if(fabs(phi_deriv(0.86)) < 1) 
        cout << "условие сходимости выполнено\n";
    else 
        cout << "НЕТ\n";
    cout << endl;

    double q_max = 0.0;
    for (double x = 0.5; x <= 1.0; x += 0.01) {
        double q = fabs(phi_deriv(x));
        if (q > q_max) q_max = q;
    }

    cout << "Максимальное |phi'(x)| на отрезке: " << q_max << endl;
    

    double x0;
    cout << "Введите начальное приближение из интервала [0.5, 1.0]: ";
    cin >> x0;
    
    double Eps = 1e-6;
    int N;
    
    cout << "================================\n\n";
    double root = ITER(x0, Eps, Delta, N);
    cout << "Корень: " << root << ", итераций: " << N << endl;
    cout << "Проверка: f(x) = " << tan(root) - 1.0/root << endl << endl;
    
    cout << "\nИССЛЕДОВАНИЕ СКОРОСТИ СХОДИМОСТИ\n";
    cout << "================================\n";

    printf("--------------------------------------------------------\n");
    printf("    Eps      | Итерации |      Корень      | Невязка\n");
    printf("--------------------------------------------------------\n");

    Eps = 0.1;
    while (Eps >= 0.00000001){
        double Root = ITER(x0, Eps, Delta, N);
        double Residual = fabs(0.860333580000 - Root);
        printf("%10.7f | %8d | %16.12f | %10.3e\n", Eps, N, Root, Residual);
        Eps /= 10.0;
    }
    
    cout << "\nИССЛЕДОВАНИЕ ЧУВСТВИТЕЛЬНОСТИ К ОШИБКАМ ОКРУГЛЕНИЯ\n";
    cout << "================================================\n";
    cout << "Delta\t\tИтерация\tКорень\t\tЗначение\n";
    cout << "----------------------------------------\n";
    
    double Delta_values[] = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11};
    int Nit = 0;
    Eps = 1e-6;
    for(double d : Delta_values) {
        double exp_root = ITER(x0, Eps, d, Nit);
        double error = Round(F(exp_root), d);
        printf("%.0e | %8d | %16.12f | %10.3e\n", 
               d, Nit, exp_root, error);
    }
    
    return 0;
}