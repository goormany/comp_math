#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double F(double x);
double BISECT(double Left, double Right, double Eps, int &N);
double Round(double X, double Delta);

double F(double x){
    return tan(x) - 1.0/x;
}

double BISECT(double Left, double Right, double Eps, int &N){
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
    if (FLeft == 0.0) return Left;
    if (FRight == 0.0) return Right;
    
    while ((Right - Left) >= E){
        X = 0.5 * (Right + Left);     /* вычисление середины отрезка */
        Y = F(X);
        if (Y == 0.0) return (X);
        if (Y * FLeft < 0.0)
            Right = X;
        else{ 
            Left = X; 
            FLeft = Y; 
        }
        N++;
    };
    
    return X;
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
    
    printf("Исследование метода бисекции для уравнения tg(x) - 1/x = 0\n");
    printf("Интервал поиска: [%g, %g]\n\n", Left, Right);
    printf("Пункт 4. Зависимость числа итераций от точности Eps:\n");
    printf("--------------------------------------------------------\n");
    printf("    Eps      | Итерации |      Корень      | Невязка\n");
    printf("--------------------------------------------------------\n");
    
    Eps = 0.1;
    while (Eps >= 0.000001){
        Root = BISECT(Left, Right, Eps, N);
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
    double TrueRoot = BISECT(Left, Right, TrueEps, N);
    
    double Delta = 0.1;
    while (Delta >= 0.000001){
        double L = Left, R = Right;
        double FL = Round(F(L), Delta);
        double FR = Round(F(R), Delta);
        int Nit = 0;
        double Xm;
        
        if (FL * FR > 0.0) { 
            printf("Delta=%g: Неверное задание интервала после округления\n", Delta);
            Delta /= 10.0;
            continue; 
        }
        
        while ((R - L) >= TrueEps){
            Xm = 0.5 * (R + L);
            double Ym = Round(F(Xm), Delta);
            
            if (Ym == 0.0) break;
            if (Ym * FL < 0.0){
                R = Xm;
                FR = Ym;
            }
            else{ 
                L = Xm; 
                FL = Ym; 
            }
            Nit++;
        }
        
        double Error = fabs(Xm - TrueRoot);
        
        printf("%10.7f | %8d | %16.12f | %10.3e\n", 
               Delta, Nit, Xm, Error);
        
        Delta /= 10.0;
    }
    
    printf("\n\nТеоретическая оценка числа итераций:\n");
    printf("N = log2((b-a)/Eps)\n");
    printf("Для Eps = 1e-6: N = log2(0.5/1e-6) = log2(500000) ≈ %d\n", (int)(log2(0.5/1e-6) + 0.5));
    
    return 0;
}