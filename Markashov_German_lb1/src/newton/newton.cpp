#include <iostream>
#include <cmath>

double delta = 0.0;

double Round (double X,double Delta){
 if (Delta<=1E-9) {std::cout << "Неверное задание точности округления" << std::endl;exit(1);}
 if (X>0.0) return (Delta*(long((X/Delta)+0.5)));
    else    return (Delta*(long((X/Delta)-0.5)));
}


double F(double x){
    double fx = (x*x - x*x*x - 1) / (4 + (x*x));
    if (::delta == 0.0){
        return fx;
    }
    return Round(fx, ::delta);
}

double F1(double x){
    double fx = (-x*x*x*x - 12 * x*x + 10 * x) / ((4 + x*x)*(4 + x*x));
    if (::delta == 0.0){
        return fx;
    }
    return Round(fx, ::delta);
}

double F2(double x){
    double fx = (8 * x*x*x - 30 * x*x - 96 * x + 40) / ((4 + x*x)*(4 + x*x)*(4 + x*x));
    if (::delta == 0.0){
        return fx;
    }
    return Round(fx, ::delta);
}

double NEWTON (double X,double Eps,int &N)
{
  extern double F1 (double);
  double Y,Y1,DX;
  N=0;
  do
  {
    Y  = F(X);
    if (Y==0.0) return (X);

    Y1 = F1(X);
    if (Y1==0.0) {puts("Производная обратилась в ноль\n");exit(1);}

    DX=Y/Y1; X=X-DX; N++;
  }
  while (fabs(DX)>Eps);
  return (X);
}

int main(){
    double left;
    double right;
    int iter;
    int x0;

    std::cout << "x0: ";
    std::cin >> x0;

    if(F(x0) * F2(x0) <= 0){
        std::cout << "Выбор начального приближения не увдоволетворяет условию" << std::endl;
        return 1;
    }

    std::cout << "Eps\tRoot\t\tIter" << std::endl;
    for(double eps = 0.1; eps >= 0.000001; eps /= 10){
        double root = NEWTON(x0, eps, iter);
        std::cout << eps << "\t" << root << "\t\t" << iter << std::endl;
    }

    std::cout << "\nDelta\tRoot\t\tIter" << std::endl;
    double eps = 0.00001;
    for(double d = 0.1; d >= 0.000001; d /= 10){
        delta = d;
        double root = NEWTON(x0, eps, iter);
        std::cout << delta << "\t" << root << "\t\t" << iter << std::endl;
    }

    return 0;
}