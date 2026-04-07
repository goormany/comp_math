#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "methods.h"

double Delta = 0.0;   // Глобальная переменная для моделирования ошибок

// Функция F(x) = e^x - arccos(sqrt(x))
double F(double x) {
    if (x < 0 || x > 1) return NAN;
    double value = exp(x) - acos(sqrt(x));
    if (Delta > 0.0)
        return Round(value, Delta);   // Моделирование ошибок
    else
        return value;   // Без округления
}

// Производная F'(x) = e^x + 1/(2*sqrt(x)*sqrt(1-x))
double F1(double x) {
    if (x <= 0 || x >= 1) return NAN;
    double value = exp(x) + 1.0/(2.0*sqrt(x)*sqrt(1.0 - x));
    if (Delta > 0.0)
        return Round(value, Delta);
    else
        return value;
}

// Вторая производная F''(x) для метода Ньютона
double F2(double x) {
    if (x <= 0 || x >= 1) return NAN;
    double sqrtX = sqrt(x);
    double sqrt1minusX = sqrt(1.0 - x);
    // F''(x) = e^x - (1-2x)/(4*x^(3/2)*(1-x)^(3/2))
    double value = exp(x) - (1.0 - 2.0*x) / (4.0 * x * sqrtX * (1.0 - x) * sqrt1minusX);
    if (Delta > 0.0)
        return Round(value, Delta);
    else
        return value;
}

// Функция Phi(x) = x - lambda*F(x) для метода простых итераций
// lambda = 0.2 (выбрано экспериментально)
double Phi(double x) {
    if (x < 0 || x > 1) return NAN;
    const double lambda = 0.2;
    double value = x - lambda * (exp(x) - acos(sqrt(x)));
    if (Delta > 0.0)
        return Round(value, Delta);
    else
        return value;
}

// Производная Phi'(x) = 1 - lambda*F'(x)
double Phi1(double x) {
    if (x <= 0 || x >= 1) return NAN;
    const double lambda = 0.2;
    double fPrime = exp(x) + 1.0/(2.0*sqrt(x)*sqrt(1.0 - x));
    double value = 1.0 - lambda * fPrime;
    if (Delta > 0.0)
        return Round(value, Delta);
    else
        return value;
}

int main() {
    double Left = 0.0;
    double Right = 1.0;

    std::ofstream csvConvergence("convergence_data.csv");
    std::ofstream csvSensitivity("sensitivity_data.csv");
    
    csvConvergence << "Method,Eps,Root,Iterations\n";
    csvSensitivity << "Method,Delta,Root,Iterations\n";

    // ============================
    // Метод бисекции
    // ============================
    std::cout << "=== Метод бисекции ===\n";

    double EpsValues[] = {0.1, 0.01, 0.001, 0.0001, 1e-5, 1e-6};
    int numEps = sizeof(EpsValues)/sizeof(EpsValues[0]);

    // --- Зависимость сходимости ---
    std::cout << "Зависимость числа итераций от Eps:\n";
    std::cout << std::left << std::setw(12) << "Eps"
              << std::setw(20) << "Корень"
              << "Число итераций\n";
    std::cout << std::string(55,'-') << "\n";

    for(int i=0;i<numEps;i++) {
        Delta = 0.0;// Без ошибок
        int iterations = 0;
        double root = BISECT(Left, Right, EpsValues[i], iterations);
        std::cout << std::left << std::setw(12) << EpsValues[i]
                  << std::setw(20) << root
                  << iterations << "\n";
        csvConvergence << "Bisect," << EpsValues[i] << "," << root << "," << iterations << "\n";
    }

    // --- Чувствительность к ошибкам ---
    double DeltaValues[] = {0.1, 0.01, 0.001, 0.0001};
    int numDelta = sizeof(DeltaValues)/sizeof(DeltaValues[0]);
    double FixedEps = 1e-6;

    std::cout << "Чувствительность метода бисекции к ошибкам (Delta):\n";
    std::cout << std::left << std::setw(12) << "Delta"
              << std::setw(20) << "Корень"
              << "Число итераций\n";
    std::cout << std::string(55,'-') << "\n";

    for(int i=0;i<numDelta;i++) {
        Delta = DeltaValues[i];
        int iterations = 0;
        double root = BISECT(Left, Right, FixedEps, iterations);
        std::cout << std::left << std::setw(12) << Delta
                  << std::setw(20) << root
                  << iterations << "\n";
        csvSensitivity << "Bisect," << DeltaValues[i] << "," << root << "," << iterations << "\n";
    }

    // ============================
    // Метод хорд
    // ============================
    std::cout << "\n=== Метод хорд ===\n";

    // --- Зависимость сходимости ---
    std::cout << "Зависимость числа итераций от Eps:\n";
    std::cout << std::left << std::setw(12) << "Eps"
              << std::setw(20) << "Корень"
              << "Число итераций\n";
    std::cout << std::string(55,'-') << "\n";

    for(int i=0;i<numEps;i++) {
        Delta = 0.0;
        int iterations = 0;
        double root = HORDA(Left, Right, EpsValues[i], iterations);
        std::cout << std::left << std::setw(12) << EpsValues[i]
                  << std::setw(20) << root
                  << iterations << "\n";
        csvConvergence << "Horda," << EpsValues[i] << "," << root << "," << iterations << "\n";
    }

    // --- Чувствительность к ошибкам ---
    std::cout << "Чувствительность метода хорд к ошибкам (Delta):\n";
    std::cout << std::left << std::setw(12) << "Delta"
              << std::setw(20) << "Корень"
              << "Число итераций\n";
    std::cout << std::string(55,'-') << "\n";

    for(int i=0;i<numDelta;i++) {
        Delta = DeltaValues[i];
        int iterations = 0;
        double root = HORDA(Left, Right, FixedEps, iterations);
        std::cout << std::left << std::setw(12) << Delta
                  << std::setw(20) << root
                  << iterations << "\n";
        csvSensitivity << "Horda," << DeltaValues[i] << "," << root << "," << iterations << "\n";
    }

    // ============================
    // Метод Ньютона
    // ============================
    std::cout << "\n=== Метод Ньютона (касательных) ===\n";

    // Условие сходимости: f(x0)*f''(x0) > 0
    // F''(x) > 0 на (0,1), корень x* ≈ 0.154
    // Берем x0 = 0.3 из области f(x) > 0, близко к корню
    double X0 = 0.3;
    
    std::cout << "Начальное приближение x0 = " << X0 << "\n";
    std::cout << "f(x0) = " << F(X0) << "\n";
    std::cout << "f'(x0) = " << F1(X0) << "\n\n";

    // --- Зависимость сходимости ---
    std::cout << "Зависимость числа итераций от Eps:\n";
    std::cout << std::left << std::setw(12) << "Eps"
              << std::setw(20) << "Корень"
              << "Число итераций\n";
    std::cout << std::string(55,'-') << "\n";

    for(int i=0;i<numEps;i++) {
        Delta = 0.0;
        int iterations = 0;
        double root = NEWTON(X0, EpsValues[i], iterations, Left, Right);
        std::cout << std::left << std::setw(12) << EpsValues[i]
                  << std::setw(20) << root
                  << iterations << "\n";
        csvConvergence << "Newton," << EpsValues[i] << "," << root << "," << iterations << "\n";
    }

    // --- Чувствительность к ошибкам ---
    std::cout << "Чувствительность метода Ньютона к ошибкам (Delta):\n";
    std::cout << std::left << std::setw(12) << "Delta"
              << std::setw(20) << "Корень"
              << "Число итераций\n";
    std::cout << std::string(55,'-') << "\n";

    for(int i=0;i<numDelta;i++) {
        Delta = DeltaValues[i];
        int iterations = 0;
        double root = NEWTON(X0, FixedEps, iterations, Left, Right);
        std::cout << std::left << std::setw(12) << Delta
                  << std::setw(20) << root
                  << iterations << "\n";
        csvSensitivity << "Newton," << DeltaValues[i] << "," << root << "," << iterations << "\n";
    }

    // ============================
    // Метод простых итераций
    // ============================
    std::cout << "\n=== Метод простых итераций ===\n";
    
    std::cout << "Проверка условия сходимости |Phi'(x)| < 1:\n";
    double maxPhiDeriv = 0.0;
    for(double x = Left; x <= Right; x += 0.01) {
        double phiDeriv = fabs(Phi1(x));
        if(phiDeriv > maxPhiDeriv) maxPhiDeriv = phiDeriv;
    }
    std::cout << "max|Phi'(x)| на [" << Left << ", " << Right << "] = " << maxPhiDeriv << "\n";
    
    if(maxPhiDeriv >= 1.0) {
        std::cout << "ВНИМАНИЕ: Условие сходимости не выполнено!\n";
    } else {
        std::cout << "Условие сходимости выполнено: q = " << maxPhiDeriv << " < 1\n";
    }
    
    double X0_iter = 0.5;
    std::cout << "\nНачальное приближение x0 = " << X0_iter << "\n";
    std::cout << "Phi(x0) = " << Phi(X0_iter) << "\n";
    std::cout << "Phi'(x0) = " << Phi1(X0_iter) << "\n\n";

    // --- Зависимость сходимости ---
    std::cout << "Зависимость числа итераций от Eps:\n";
    std::cout << std::left << std::setw(12) << "Eps"
              << std::setw(20) << "Корень"
              << "Число итераций\n";
    std::cout << std::string(55,'-') << "\n";

    for(int i=0;i<numEps;i++) {
        Delta = 0.0;
        int iterations = 0;
        double root = ITER(X0_iter, EpsValues[i], iterations);
        std::cout << std::left << std::setw(12) << EpsValues[i]
                  << std::setw(20) << root
                  << iterations << "\n";
        csvConvergence << "Iter," << EpsValues[i] << "," << root << "," << iterations << "\n";
    }

    // --- Чувствительность к ошибкам ---
    std::cout << "Чувствительность метода простых итераций к ошибкам (Delta):\n";
    std::cout << std::left << std::setw(12) << "Delta"
              << std::setw(20) << "Корень"
              << "Число итераций\n";
    std::cout << std::string(55,'-') << "\n";

    for(int i=0;i<numDelta;i++) {
        Delta = DeltaValues[i];
        int iterations = 0;
        double root = ITER(X0_iter, FixedEps, iterations);
        std::cout << std::left << std::setw(12) << Delta
                  << std::setw(20) << root
                  << iterations << "\n";
        csvSensitivity << "Iter," << DeltaValues[i] << "," << root << "," << iterations << "\n";
    }

    csvConvergence.close();
    csvSensitivity.close();
    
    std::cout << "\nДанные сохранены в convergence_data.csv и sensitivity_data.csv\n";

    Delta = 0.0;
    return 0;
}
