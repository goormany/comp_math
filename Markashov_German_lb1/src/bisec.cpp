#include <iostream>
#include <cmath>
#include <iomanip>

double delta = 0.0;
double Round (double X,double Delta);

double F(double x){
    double fx = (x*x - x*x*x - 1) / (4 + (x*x));
    if (::delta == 0.0){
        return fx;
    }
    return Round(fx, ::delta);
}

double f_p(double x){
    return (-(x*x*x*x) - 12 * (x*x) + 10 * x) / ((4 + (x*x)) * (4 + (x*x)));
}


double BISECT(double Left,double Right,double Eps,int &N)
{
  double E = fabs(Eps)*2.0;
  double FLeft = F(Left);
  double FRight = F(Right);
  double X = (Left+Right)/2.0;
  double Y;

  if (FLeft*FRight>0.0) {puts("Неверное задание интервала\n");exit(1);}
  if (Eps<=0.0) {puts("Неверное задание точности\n");exit(1);}

  N=0;
  if (FLeft==0.0)  return Left;
  if (FRight==0.0) return Right;

  while ((Right-Left)>=E)
  {
	X = 0.5*(Right + Left);     /* вычисление середины отрезка       */
	Y = F(X);
	if (Y == 0.0) return (X);
        if (Y*FLeft < 0.0)
           Right=X;
        else
         { Left=X; FLeft=Y; }
        N++;
  };
  return(X);
}

double Round (double X,double Delta)
{
 if (Delta<=1E-9) {std::cout << "Неверное задание точности округления" << std::endl;exit(1);}
 if (X>0.0) return (Delta*(long((X/Delta)+0.5)));
    else    return (Delta*(long((X/Delta)-0.5)));
}

int main(){
    double left = -1;
    double right = 1;
    int iter;

    std::cout << "Зависимость числа итераций от точности" << std::endl;
    std::cout << "Eps\tRoot\t\tIter" << std::endl;

    for(double eps = 0.1; eps >= 0.000001; eps /= 10){
        double root = BISECT(left, right, eps, iter);
        std::cout << eps << "\t" << root << "\t\t" << iter << std::endl;
    }

    std::cout << "Исследование ошибок" << std::endl;
    std::cout << "Delta\tRoot\t\tIter" << std::endl;
    iter = 0;
    double eps = 0.00001;
    for(double d = 0.1; d >= 0.000001; d /= 10){
        delta = d;
        double root = BISECT(left, right, eps, iter);
        std::cout << delta << "\t" << root << "\t\t" << iter << std::endl;
    }

    return 0;
}