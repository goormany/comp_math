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
        return x + 2.173 * fx;
    }

    fx = x + 2.173 * fx;

    return Round(fx, ::delta);
}


double ITER(double X0,double Eps,int &N)
{
  if (Eps<=0.0) {puts("Неверное задание точности\n");exit (1);}
  double X1=F(X0);
  double X2=F(X1);
  N = 2;
  while( (X1 - X2)*(X1 - X2) > fabs((2*X1-X0-X2)*Eps) )
  {
	X0 = X1;
	X1 = X2;
	X2 = F(X1);
        N++;
  }
   return(X2);
}

int main(){
    double left;
    double right;
    int iter;
    double x0;

    std::cout << "x0: ";
    std::cin >> x0;


    std::cout << "Eps\tRoot\t\tIter" << std::endl;
    for(double eps = 0.1; eps >= 0.000001; eps /= 10){
        double root = ITER(x0, eps, iter);
        std::cout << eps << "\t" << root << "\t\t" << iter << std::endl;
    }

    std::cout << "\nDelta\tRoot\t\tIter" << std::endl;
    double eps = 0.00001;
    for(double d = 0.001; d >= 0.000001; d /= 10){
        delta = d;
        double root = ITER(x0, eps, iter);
        std::cout << delta << "\t" << root << "\t\t" << iter << std::endl;
    }

    return 0;
}