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

double f(double x) {
    return 2*x*x - x*x*x*x - 1 - log(x);
}

double df(double x) {
    return 4*x - 4*x*x*x - 1/x;
}

double Phi(double x) {
    double lambda_opt = -0.609;
    return x - lambda_opt * f(x);
}

double dphi(double x) {
    double lambda_opt = -0.609;
    return 1.0 - lambda_opt * df(x);
}

double F_round(double x, double delta) {
    return Round(Phi(x), delta);
}

double ITER(double X0,double Eps,int &N)
{
  if (Eps<=0.0) {puts("Неверное задание точности\n");exit (1);}
  double X1=Phi(X0);
  double X2=Phi(X1);
  N = 2;
  while( (X1 - X2)*(X1 - X2) > fabs((2*X1-X0-X2)*Eps) )
  {
	X0 = X1;
	X1 = X2;
	X2 = Phi(X1);
    N++;
  }
   return(X2);
}

int main() {
    double x0 = 1.15;


    double eps = 1e-6;
    double Delta = 0.000001;
    int iterations;
    double root;

    // cout << "Введите начальное значение x0: ";
    // cin >> x0;

    // cout << "Введите точность ответа: ";
    // cin >> eps;

    cout << "МЕТОД ПРОСТЫХ ИТЕРАЦИЙ\n";
    cout << "Ответ = " << ITER(x0, eps, iterations) << endl;
    cout << "Количество итераций = " << iterations << endl << endl;


    cout << "ЧАСТЬ 1: СКОРОСТЬ СХОДИМОСТИ\n";
    cout << "Eps\t\tКорень\t\tИтераций\t|F(x)|\n";
    cout << "------------------------------------------------\n";

    eps = 0.1;
    for (int i = 0; i < 10; i++) {
        root = ITER(x0, eps, iterations);
        cout << eps << "\t\t" << root << "\t\t" << iterations << "\t\t" << fabs(f(root)) << endl;
        eps /= 10;
    }



    cout << "\nЧАСТЬ 2: ИССЛЕДОВАНИЕ ОБУСЛОВЛЕННОСТИ\n";
    cout << "------------------------------------------------\n";
    double delta = eps;

    double exact_root = ITER(x0, 1e-10, iterations);
    double nu = 1.0 / fabs(dphi(exact_root));
    
    cout << "Точный корень: " << exact_root << endl;
    cout << "Число обусловленности ν_Δ = " << nu << endl;
    cout << "ν_Δ * Delta = " << nu * delta << endl;

    #define Phi(x) F_round(x, delta)
    
    double round_root = ITER(x0, 1e-9, iterations);
    double error = fabs(exact_root - round_root);
    
    cout << "Корень с округлением: " << round_root << endl;
    cout << "Разница |x_точн - x_округл| = " << error << endl;
    
    return 0;
}
