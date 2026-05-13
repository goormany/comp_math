#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace std;

// Исходная функция
double f(double x) {
    return atan(x) - 1.0 / x;
}

// Первая производная
double df(double x) {
    return 1.0 / (1.0 + x * x) + 1.0 / (x * x);
}

// Барицентрическая форма Лагранжа
vector<double> barycentricWeights(const vector<double>& x) {
    const int n = (int)x.size();
    vector<double> w(n, 1.0);

    for (int i = 0; i < n; ++i) {
        long double prod = 1.0L;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            prod *= (x[i] - x[j]);
        }
        w[i] = 1.0 / (double)prod;
    }

    // Нормировка весов
    double mx = 0.0;
    for (double wi : w) mx = max(mx, fabs(wi));
    if (mx > 0) {
        for (double& wi : w) wi /= mx;
    }

    return w;
}

double lagrangeBarycentric(const vector<double>& x,
                           const vector<double>& y,
                           const vector<double>& w,
                           double xq) {
    const double eps = 1e-14;
    const int n = (int)x.size();

    for (int i = 0; i < n; ++i) {
        if (fabs(xq - x[i]) < eps) return y[i];
    }

    long double num = 0.0L, den = 0.0L;
    for (int i = 0; i < n; ++i) {
        long double t = (long double)w[i] / (xq - x[i]);
        num += t * y[i];
        den += t;
    }
    return (double)(num / den);
}

// Метод прогонки
// Решает трёхдиагональную систему:
// a[i] * x[i-1] + b[i] * x[i] + c[i] * x[i+1] = d[i]
vector<double> solveTridiagonal(const vector<double>& a,
                                const vector<double>& b,
                                const vector<double>& c,
                                const vector<double>& d) {
    const int n = (int)b.size();
    vector<double> cp(n, 0.0), dp(n, 0.0), x(n, 0.0);

    cp[0] = (n > 1 ? c[0] / b[0] : 0.0);
    dp[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double denom = b[i] - a[i] * cp[i - 1];
        if (i < n - 1) cp[i] = c[i] / denom;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }

    x[n - 1] = dp[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    return x;
}

// Натуральный кубический сплайн
// M[0] = M[n-1] = 0
// M[i] - вторые производные в узлах
struct CubicSpline {
    vector<double> x;   // узлы
    vector<double> y;   // значения в узлах
    vector<double> M;   // моменты - вторые производные

    double eval(double xq) const {
        const int n = (int)x.size();
        if (xq <= x.front()) return y.front();
        if (xq >= x.back())  return y.back();

        auto it = upper_bound(x.begin(), x.end(), xq);
        int i = (int)(it - x.begin()) - 1;

        double h = x[i + 1] - x[i];
        double A = (x[i + 1] - xq) / h;
        double B = (xq - x[i]) / h;

        double S =
            M[i]     * pow(x[i + 1] - xq, 3) / (6.0 * h) +
            M[i + 1] * pow(xq - x[i],     3) / (6.0 * h) +
            (y[i]     - M[i]     * h * h / 6.0) * A +
            (y[i + 1] - M[i + 1] * h * h / 6.0) * B;

        return S;
    }
};

// Закреплённый кубический сплайн
// Граничные условия: S'(a) = f'(a), S'(b) = f'(b)
CubicSpline buildClampedCubicSpline(const vector<double>& x,
                                    const vector<double>& y) {
    const int n = (int)x.size();
    CubicSpline spline;
    spline.x = x;
    spline.y = y;
    spline.M.assign(n, 0.0);

    if (n <= 2) return spline; // сплайн-интерполяция не определена

    const int N = n;
    vector<double> a(N, 0.0), b(N, 0.0), c(N, 0.0), d(N, 0.0);

    // Левая граница (i = 0)
    double h1 = x[1] - x[0];
    b[0] = h1 / 3.0;
    c[0] = h1 / 6.0;
    d[0] = (y[1] - y[0]) / h1 - df(x[0]);

    // Внутренние узлы (i = 1 ... n-2)
    for (int i = 1; i <= n - 2; ++i) {
        double hPrev = x[i] - x[i - 1];
        double hNext = x[i + 1] - x[i];

        a[i] = hPrev / 6.0;
        b[i] = (hPrev + hNext) / 3.0;
        c[i] = hNext / 6.0;
        d[i] = (y[i + 1] - y[i]) / hNext - (y[i] - y[i - 1]) / hPrev;
    }

    // Правая граница (i = n-1)
    double hLast = x[n - 1] - x[n - 2];
    a[n - 1] = hLast / 6.0;
    b[n - 1] = hLast / 3.0;
    d[n - 1] = df(x[n - 1]) - (y[n - 1] - y[n - 2]) / hLast;

    // Решаем трёхдиагональную систему
    vector<double> M_all = solveTridiagonal(a, b, c, d);

    // Копируем все моменты
    for (int i = 0; i < n; ++i) {
        spline.M[i] = M_all[i];
    }

    return spline;
}

// Таблица для дальнейшего анализа
struct Sample {
    double x;
    double fx;
    double px;
    double sx;
    double errP;
    double errS;
};

int main() {
    const double A = 0.1;
    const double B = 5.0;

    int n;
    cout << "Введите число узлов n: ";
    cin >> n;

    if (n < 3) {
        cerr << "Ошибка: n должно быть >= 3\n";
        return 1;
    }

    // Равномерные узлы
    vector<double> x(n), y(n);
    double h = (B - A) / (n - 1);

    for (int i = 0; i < n; ++i) {
        x[i] = A + i * h;
        y[i] = f(x[i]);
    }

    ofstream nodesOut("nodes.csv");
    nodesOut << "x,y\n";
    for (int i = 0; i < n; ++i) {
        nodesOut << setprecision(15) << x[i] << "," << y[i] << "\n";
    }
    nodesOut.close();

    // Весы для барицентрической формулы Лагранжа
    vector<double> w = barycentricWeights(x);

    // Сплайн
    CubicSpline spline = buildClampedCubicSpline(x, y);

    // Печать узлов
    cout << fixed << setprecision(10);
    cout << "\nУзлы интерполяции:\n";
    cout << "i\t x_i\t\t f(x_i)\n";
    for (int i = 0; i < n; ++i) {
        cout << i << "\t " << x[i] << "\t " << y[i] << "\n";
    }

    // Сетка для анализа
    const int Nfine = 1000;
    vector<Sample> samples;
    samples.reserve(Nfine + 1);

    double maxErrP = 0.0;
    double maxErrS = 0.0;

    for (int i = 0; i <= Nfine; ++i) {
        double t = (double)i / Nfine;
        double xx = A + (B - A) * t;

        double fx = f(xx);
        double px = lagrangeBarycentric(x, y, w, xx);
        double sx = spline.eval(xx);

        double errP = fabs(px - fx);
        double errS = fabs(sx - fx);

        maxErrP = max(maxErrP, errP);
        maxErrS = max(maxErrS, errS);

        samples.push_back({xx, fx, px, sx, errP, errS});
    }

    cout << "\nМаксимальная ошибка на мелкой сетке:\n";
    cout << "Lagrange = " << maxErrP << "\n";
    cout << "Spline   = " << maxErrS << "\n";

    // Сохранение данных для построения графиков в Python
    ofstream out("interpolation_data.csv");
    out << "x,f(x),lagrange,spline,err_lagrange,err_spline\n";
    for (const auto& s : samples) {
        out << setprecision(15)
            << s.x << ","
            << s.fx << ","
            << s.px << ","
            << s.sx << ","
            << s.errP << ","
            << s.errS << "\n";
    }
    out.close();

    cout << "\nДанные сохранены в файл interpolation_data.csv\n";

    return 0;
}