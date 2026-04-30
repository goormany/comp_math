#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "Task.h"


using namespace std;

double Delta = 1e-4;
bool Perturb = false;

double F(double x) {
  double val = sin(x * x) - 6.0 * x + 1.0;

  if (Perturb && Delta > 0) {
    val = Round(val, Delta);
  }
  return val;
}

double F1(double x) {
  double val = 2.0 * x * cos(x * x) - 6.0;
  
  if (Perturb && Delta > 0) {
    val = Round(val, Delta);
  }
  return val;
}

int main() {

  cout << "Lab 3-4: Bisect's method" << endl;
  cout << "f(x) = sin(x^2) - 6x + 1" << endl;

  double left, right;
  int N = 0;

  std::cout << "left: ";
  std::cin >> left;


  std::cout << "right: ";
  std::cin >> right;

  if (F(left) * F(right) > 0) {
      cout << "На интервале [" << left << ", " << right 
           << "] нет гарантии наличия корня! f(left)*f(right) > 0." << endl;
      return 1;
  }

  double x0 = (left + right) / 2.0;

  cout << "\nBisect - dependence of N on Eps:" << endl;
  cout << "Eps; N; X; f(X)" << endl;

  for (int i = 0; i <= 6; i++) {
    Perturb = false;
    double Eps = pow(10, -i);
    double x = BISECT(left, right, Eps, N);
    cout << scientific << Eps << "; " << N << "; " << fixed << setprecision(8) << x << "; " << F(x) << endl;
  }

  cout << "\n" << endl;
  cout << "Delta; X; N" << endl;
  double Eps = 0.00001;
  for (int i = 0; i <= 8; i++) {
    double delta = pow(10, -i);
    N = 0; 
    Perturb = true;
    Delta = delta; 
    double x = BISECT(left, right, Delta, N);
    cout << scientific << Delta << "; " << fixed << setprecision(8) << x << "; " << N << endl;
  
  }
  Perturb = false;
  Delta = 1e-4;
  
cout << "\nHORD - dependence of N on Eps:" << endl;
cout << "Eps; N; X; f(X)" << endl;

if (F(left) * F(right) > 0) {
    std::cout << "На интервале [" << left << ", " << right 
              << "] нет гарантии наличия корня!" << std::endl;
    std::cout << "F(left) = " << F(left) << ", F(right) = " << F(right) << std::endl;
    return 1;
}

for (int i = 0; i <= 6; i++) {
    double eps = pow(10, -i);
    N = 0;
    Perturb = false;
    double x = HORD(left, right, Eps, N);
    cout << scientific << Delta << "; " << fixed << setprecision(8) << x << "; " << N << endl;
  }

cout << "\n" << endl;

double eps = 0.00001;
cout << "Delta; X; N" << endl;

for (int i = 0; i <= 8; i++) {
    double delta = pow(10, -i);
    N = 0; 
    Perturb = true;
    Delta = delta; 
    double x = HORD(left, right, eps, N);
    cout << scientific << delta << "; " << fixed << setprecision(8) << x << "; " << N << endl;
}

cout << "\nDone! Excel." << endl;

return 0;
}
