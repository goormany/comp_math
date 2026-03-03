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

double HORDA(double Left,double Right,double Eps,int &N){
  double FLeft = F(Left);
  double FRight = F(Right);
  double X,Y;

  if (FLeft*FRight>0.0) {puts("Неверное задание интервала\n");exit(1);}
  if (Eps<=0.0) {puts("Неверное задание точности\n");exit(1);}

  N=0;
  if (FLeft==0.0)  return Left;
  if (FRight==0.0) return Right;

  do
  {
	X = Left-(Right-Left)*FLeft/(FRight-FLeft);
	Y = F(X);
	if (Y == 0.0) return (X);
        if (Y*FLeft < 0.0)
         { Right=X; FRight=Y; }
        else
         { Left=X; FLeft=Y; }
        N++;
  }
  while ( fabs(Y) >= Eps );

  return(X);

}


int main(){
    double left;
    double right;
    int iter;

    std::cout << "left: ";
    std::cin >> left;

    std::cout << "right: ";
    std::cin >> right;

    if (F(left) * F(right) > 0) {
        std::cout << "На интервале [" << left << ", " << right 
                  << "] нет гарантии наличия корня!" << std::endl;
        std::cout << "F(left) = " << F(left) << ", F(right) = " << F(right) << std::endl;
        return 1;
    }

    std::cout << "Eps\tRoot\t\tIter" << std::endl;
    for(double eps = 0.1; eps >= 0.000001; eps /= 10){
        double root = HORDA(left, right, eps, iter);
        std::cout << eps << "\t" << root << "\t\t" << iter << std::endl;
    }

    std::cout << "Исследование ошибок" << std::endl;
    std::cout << "Delta\tRoot\t\tIter" << std::endl;
    double eps = 0.00001;
    for(double d = 0.1; d >= 0.000001; d /= 10){
        delta = d;
        double root = HORDA(left, right, eps, iter);
        std::cout << delta << "\t" << root << "\t\t" << iter << std::endl;

    }

    return 0;
}