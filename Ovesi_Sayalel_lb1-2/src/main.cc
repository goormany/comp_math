
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

double delta_global = 0.0;

double Round(double x, double delta) {
  if (delta <= 1e-9) {
    std::cerr << "Round: delta too small\n";
    std::exit(1);
  }
  if (x > 0.0)
    return delta * (long)((x / delta) + 0.5);
  else
    return delta * (long)((x / delta) - 0.5);
}

// f(x) = 2^(x^2) - 1/x
double F(double x) {
  if (x == 0.0) {
    std::cerr << "F: division by zero\n";
    std::exit(1);
  }
  double val = std::pow(2.0, x * x) - 1.0 / x;
  if (delta_global > 0.0)
    val = Round(val, delta_global);
  return val;
}

// f'(x) = 2^(x^2) * ln(2) * 2x + 1/x^2
double F1(double x) {
  if (x == 0.0) {
    std::cerr << "F1: division by zero\n";
    std::exit(1);
  }
  double val = std::pow(2.0, x * x) * std::log(2.0) * 2.0 * x + 1.0 / (x * x);
  if (delta_global > 0.0)
    val = Round(val, delta_global);
  return val;
}

// phi(x)
//   f(x) = 0  =>  2^(x^2) = 1/x  =>  x = 2^(-x^2)
//   phi(x) = 2^(-x^2)
//
//   |phi'(x)| = 2*|x|*ln(2)*2^(-x^2) ~= 0.69 < 1  near root => converges
double PHI(double x) {
  double val = std::pow(2.0, -x * x);
  if (delta_global > 0.0)
    val = Round(val, delta_global);
  return val;
}

//   x_{n+1} = x_n - f(x_n) / f'(x_n)
//   Stop when |dx| < eps
double NEWTON(double x0, double eps, int &n) {
  if (eps <= 0.0) {
    std::cerr << "NEWTON: eps must be positive\n";
    std::exit(1);
  }

  double x = x0, dx;
  n = 0;

  do {
    double y = F(x);
    if (y == 0.0)
      return x;

    double y1 = F1(x);
    if (y1 == 0.0) {
      std::cerr << "NEWTON: derivative is zero at x = " << x << "\n";
      std::exit(1);
    }

    dx = y / y1;
    x -= dx;
    ++n;
  } while (std::fabs(dx) > eps);

  return x;
}

//   x_{n+1} = phi(x_n)
//   Steffensen stop criterion
//     (x1-x2)^2  <=  |(2*x1 - x0 - x2)| * eps
double ITER(double x0, double eps, int &n, int maxn = 500) {
  if (eps <= 0.0) {
    std::cerr << "ITER: eps must be positive\n";
    std::exit(1);
  }

  double x1 = PHI(x0);
  double x2 = PHI(x1);
  n = 2;

  while ((x1 - x2) * (x1 - x2) > std::fabs((2.0 * x1 - x0 - x2) * eps) &&
         n < maxn) {
    x0 = x1;
    x1 = x2;
    x2 = PHI(x1);
    ++n;
  }

  return x2;
}

int main() {
  std::cout << std::fixed << std::setprecision(10);

  std::cout << "f(x)  = 2^(x^2) - 1/x\n";
  std::cout << "f'(x) = 2^(x^2)*ln(2)*2x + 1/x^2\n";
  std::cout << "phi(x)= 2^(-x^2)\n";
  std::cout << "Search interval: [0.70, 0.71]\n\n";

  const double X0 = 0.70; // f(x0)*f''(x0) > 0  (verified)
  const double EPS = 1e-6;

  // 3.4  Newton's method
  std::cout << "\033[4m\033[1mNewton's method (baseline, delta=0) \033[0m\n";
  {
    int n = 0;
    double root = NEWTON(X0, EPS, n);
    std::cout << "  Root       : " << root << "\n";
    std::cout << "  f(root)    : " << F(root) << "\n";
    std::cout << "  Iterations : " << n << "\n";
  }

  std::cout << "\n\n  \033[1m\033[4mConvergence (iterations vs eps):\033[0m\n";
  std::cout << "  \033[35m\033[4m" << std::setw(14) << "eps" << std::setw(18)
            << "root" << std::setw(10) << "iters" << "\033[0m\n";
  for (double e : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6}) {
    int ni = 0;
    delta_global = 0.0;
    double r = NEWTON(X0, e, ni);
    std::cout << "  " << std::setw(14) << e << std::setw(18) << r
              << std::setw(10) << ni << "\n";
  }

  std::cout
      << "\n  \033[1m\033[4mSensitivity to input errors (delta):\033[0m\n";
  std::cout << " \033[35m\033[4m" << std::setw(14) << "delta" << std::setw(18)
            << "root" << std::setw(10) << "iters" << "\033[0m\n";
  for (double d : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6}) {
    int ni = 0;
    delta_global = d;
    double r = NEWTON(X0, EPS, ni);
    std::cout << "  " << std::setw(14) << d << std::setw(18) << r
              << std::setw(10) << ni << "\n";
  }
  delta_global = 0.0;

  // 3.5  Simple iterations
  std::cout
      << "\n\033[4m\033[1m 3.5  Simple iterations  phi(x)=2^(-x^2)\033[0m  \n";
  {
    int n = 0;
    double root = ITER(X0, EPS, n);
    std::cout << "  Root       : " << root << "\n";
    std::cout << "  f(root)    : " << F(root) << "\n";
    std::cout << "  Iterations : " << n << "\n";
  }

  std::cout << "\n  \033[4m\033[1mConvergence (iterations vs eps):\033[0m\n";
  std::cout << "  \033[35m\033[4m" << std::setw(14) << "eps" << std::setw(18)
            << "root" << std::setw(10) << "iters" << "\033[0m\n";
  for (double e : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6}) {
    int ni = 0;
    delta_global = 0.0;
    double r = ITER(X0, e, ni);
    std::cout << "  " << std::setw(14) << e << std::setw(18) << r
              << std::setw(10) << ni << "\n";
  }

  std::cout
      << "\n  \033[4m\033[1mSensitivity to input errors (delta):\n\033[0m";
  std::cout << "  \033[35m\033[4m" << std::setw(14) << "delta" << std::setw(18)
            << "root" << std::setw(10) << "iters" << "\033[0m\n";
  for (double d : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6}) {
    int ni = 0;
    delta_global = d;
    double r = ITER(X0, EPS, ni);
    std::cout << "  " << std::setw(14) << d << std::setw(18) << r
              << std::setw(10) << ni << "\n";
  }
  delta_global = 0.0;

  return 0;
}
