#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

// Подынтегральная функция: f(x) = cos(x^3)
double f(double x) {
    return cos(x * x * x);
}

// Составная формула прямоугольников (средних)
double rectangle(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double x_mid = a + (i + 0.5) * h;
        sum += f(x_mid);
    }
    return h * sum;
}

// Составная формула трапеций
double trapezoid(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; ++i) {
        sum += f(a + i * h);
    }
    return h * sum;
}

// Составная формула Симпсона
double simpson(double a, double b, int n) {
    // Для метода Симпсона n должно быть чётным
    if (n % 2 != 0) {
        n++; // если нечётное, увеличиваем на 1
    }
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        if (i % 2 == 0) {
            sum += 2.0 * f(x);  // чётные узлы
        } else {
            sum += 4.0 * f(x);  // нечётные узлы
        }
    }
    return h * sum / 3.0;
}

// Функция автоматического подбора шага с правилом Рунге
void runge_method(double a, double b, double eps, const string& name,
                  double (*method)(double, double, int), int order) {
    // Начальное число разбиений
    int n = 2;
    
    // Знаменатель для правила Рунге: 2^p - 1
    double denominator = pow(2.0, order) - 1.0;
    
    double I_h, I_h2, error;
    const int max_iter = 30;  // максимальное число удвоений
    
    for (int iter = 0; iter < max_iter; ++iter) {
        I_h = method(a, b, n);      // приближение с шагом h
        I_h2 = method(a, b, 2 * n); // приближение с шагом h/2
        
        // Оценка погрешности по правилу Рунге
        error = fabs(I_h2 - I_h) / denominator;
        
        // Если достигнута требуемая точность, выходим
        if (error <= eps) {
            break;
        }
        
        // Иначе удваиваем число разбиений
        n *= 2;
    }
    
    // Уточнённое значение интеграла (экстраполяция Ричардсона)
    double I_improved = I_h2 + (I_h2 - I_h) / denominator;
    
    // Финальное число разбиений (последнее, с которым считали I_h2)
    int n_final = 2 * n;
    double h_final = (b - a) / n_final;
    
    // Вывод результатов
    cout << fixed << setprecision(12);
    cout << name << ":\n";
    cout << "  Количество разбиений n = " << n_final << "\n";
    cout << "  Шаг интегрирования h   = " << h_final << "\n";
    cout << "  Уточнённое значение    = " << I_improved << "\n";
    cout << "  Оценка погрешности Δ   = " << scientific << setprecision(6) << error << "\n\n";
}

int main() {
    system("chcp 65001 > nul");
    
    double a = 0.0;      // нижний предел интегрирования
    double b = 1.0;      // верхний предел интегрирования
    double eps;
    
    cout << "============================================================\n";
    cout << "  Численное интегрирование функции f(x) = cos(x^3)\n";
    cout << "  на отрезке [0, 1]\n";
    cout << "============================================================\n\n";
    
    cout << "Введите требуемую точность epsilon (например, 0.01, 0.001, 0.0001): ";
    cin >> eps;
    
    cout << "\nТочность epsilon = " << eps << "\n\n";
    
    // Запуск методов с автоматическим подбором шага
    runge_method(a, b, eps, "Метод прямоугольников (средних)", rectangle, 2);
    runge_method(a, b, eps, "Метод трапеций", trapezoid, 2);
    runge_method(a, b, eps, "Метод Симпсона", simpson, 4);
    
    return 0;
}