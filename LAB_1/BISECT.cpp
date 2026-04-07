#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double f(double x) {
    return 2*x*x-x*x*x*x-1-log(x);
}

double df(double x) {
    return 4*x-4*x*x*x-1/x;
}

double Round(double x, double delta) {
    return round(x / delta) * delta;
}


double BISECT(double a, double b, double eps, int& k) {
    double c;
    k = 0;
    
    if (f(a) * f(b) >= 0) {
        cout << "Ошибка: на концах отрезка одинаковые знаки!" << endl;
        return 0;
    }
    
    while ((b - a) >= 2 * eps) {
        c = (a + b) / 2;
        
        if (f(c) == 0) {
            k++;
            return c;
        }
        
        if (f(a) * f(c) < 0) {
            b = c;
        } else {
            a = c;
        }
        k++;
    }
    
    return (a + b) / 2;
}

vector<double> BISECT_sequence(double a, double b, int N) {
    vector<double> seq;
    double left = a, right = b;
    
    for (int i = 0; i < N; i++) {
        double x = (left + right) / 2;
        seq.push_back(x);
        
        if (f(x) == 0) return seq;
        
        if (f(left) * f(x) < 0) {
            right = x;
        } else {
            left = x;
        }
    }
    return seq;
}

vector<double> BISECT_sequence_noise(double a, double b, int N, double delta) {
    vector<double> seq;
    double left = a, right = b;
    
    for (int i = 0; i < N; i++) {
        double x = (left + right) / 2;
        seq.push_back(x);
        
        double fx = Round(f(x), delta);
        double fleft = Round(f(left), delta);
        
        if (fx == 0) return seq;
        
        if (fleft * fx < 0) {
            right = x;
        } else {
            left = x;
        }
    }
    return seq;
}


int main() {
    double a = 0.99, b = 1.15;
    double eps;
    double delta = 0.00001;
    int k;
    double root;

    // cout << "Введите левую границу a: ";
    // cin >> a;

    // cout << "Введите правую границу b: ";
    // cin >> b;

    // cout << "Введите точность округления delta: ";
    // cin >> delta;

    cout << "МЕТОД БИСЕКЦИИ\n";

    
    cout << "Точность\tКорень\t\tКол-во итераций\t\tΔx\t\tmax_Δx\t\tΔx < max_Δx" << endl;
    cout << "----------------------------------------------------------------------------------------------------" << endl;
    
    eps = 0.1;
    for (int i = 0; i < 10; i++) {
        root = BISECT(a, b, eps, k);
        double ν_Δ = 1/abs(df(root));
        double max_Δx = eps / delta;
        bool check = ν_Δ < max_Δx;
        cout << eps << "\t\t" << root << "\t\t" << k << "\t\t\t" << ν_Δ << "\t" << max_Δx << "\t\t" << (check ? "ДА" : "НЕТ") << endl;
        eps /= 10;
    }



    // cout << "\n\nВТОРАЯ ЧАСТЬ: Исследование обусловленности" << endl;
    // cout << "i\tТочный x_i\tС шумом x_i*\tΔ_i = |x_N - x_i*|\tΔ_i ≤ ν·ε?" << endl;
    // cout << "----------------------------------------------------------------------------------------------------" << endl;
    
    
    // vector<double> exact_seq = BISECT_sequence(a, b, N);
    // vector<double> noise_seq = BISECT_sequence_noise(a, b, N, delta);
    // double x_n = exact_seq.back();
    // double ν_Δ = 1/abs(df(x_n));
    // double epsilon = delta * ν_Δ;

    // for (int i = 0; i < N; i++) {
    //     double delta_i = abs(x_n - noise_seq[i]);
        
    //     cout << i+1 << "\t" << exact_seq[i] << "\t" << noise_seq[i] << "\t\t" << delta_i << "\t\t" << (delta_i <= epsilon ? "ДА" : "НЕТ") << endl;

    // }
    
       
    return 0;
}