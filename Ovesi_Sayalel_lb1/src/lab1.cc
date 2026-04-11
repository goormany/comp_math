#include <cfloat>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace std;
const int T = 4;
const int M = 2;
const double EPSILON = 0.001;
double f(double x) { return pow(2.0, x * x) - 1.0f / x; }
/*void distribution() {
  cout << "\033[37m\e[4mFloating Point Distribution\e[0m\n" << endl;
  cout << "Float characteristics:" << endl;
  cout << "\tSize: " << sizeof(float) << " bytes\n";
  cout << "\tMantissa bits: " << FLT_MANT_DIG << endl;
  cout << "\tExponent bits: " << (sizeof(float) * 8 - FLT_MANT_DIG) << endl;
  cout << "\tMin exponent: " << FLT_MIN_EXP << endl;
  cout << "\tMax exponent: " << FLT_MAX_EXP << endl;
  cout << "\tMin positive: " << scientific << FLT_MIN << endl;
  cout << "\tMax: " << FLT_MAX << endl;

  cout << "Double characteristics:" << endl;
  cout << "\tSize: " << sizeof(double) << " bytes\n";
  cout << "\tMantissa bits: " << DBL_MANT_DIG << endl;
  // HOW DOES THIS WORK?
  cout << "\tExponent bits: " << (sizeof(double) * 8 - DBL_MANT_DIG) << endl;
  cout << "\tMin exponent: " << DBL_MIN_EXP << endl;
  cout << "\tMax exponent: " << DBL_MAX_EXP << endl;
  cout << "\tMin positive: " << scientific << DBL_MIN << endl;
  cout << "\tMax: " << DBL_MAX << endl << endl;

  cout << "Non uniform distribution" << endl;
  float f = 1.0f;
  for (int i = 0; i < 5; i++) {
    float next = nextafterf(f, FLT_MAX);
    cout << "Around" << f << ": spacing = " << scientific << (next - f)
         << ", relative = " << (next - f) / f << endl;
  }
}*/
void distribution() {
  cout << " \n\033[37m\e[4m\e[1m1. DISTRIBUTION OF FLOATING POINT NUMS\e[0m\n"
       << endl;
  cout << "Float params:" << std::endl;
  cout << "\tMantissa digits (T): " << T << endl;
  cout << "\tExponent digits (M): " << M << endl;
  cout << "\tBase (b): 2\n" << endl;

  // range
  int L = -(1 << (M - 1));        //-2^(m-1)
  int M_max = (1 << (M - 1)) - 1; // 2^(m-1)-1
  double t_eps = pow(2.0, -T);

  cout << "\033[37m\e[4mCharacteristics: \033[0m" << endl;
  cout << "\tExponent range: [" << L << "," << M_max << "]" << endl;
  cout << "\tMachine epsilon: " << scientific << t_eps << endl;

  // smallest positive (d1=1,o=0,n=L)
  double min_pos = pow(2.0, L) * 0.5;
  cout << "\tSmallest positive: " << scientific << min_pos << endl;
  double f_max = 0.0;
  for (int i = 0; i <= T; i++) {
    f_max += 1.0 / pow(2.0, i);
  }
  double max_num = f_max * pow(2.0, M_max);
  cout << "\tLargest number: " << scientific << max_num << endl << endl;

  cout << "\e[1m\033[37m\033[4mNon uniform distribution:\33[0m" << endl;
  cout << "(Numbers equally spaced in each exponent interval, but spacing "
          "doubles with each increase in exponent)"
       << endl
       << endl;
  double exponents[] = {static_cast<double>(L), 0.1,
                        static_cast<double>(M_max - 1)};
  const char *desc[] = {"smallest", "around 1", "around 2", "largest"};
  for (int i = 0; i < 4; i++) {
    double n = exponents[i];
    double spacing = pow(2.0, -T) * pow(2.0, n); // eps * 2^n
    double value = pow(2.0, n) * 0.5;            // approx val with d1=1
    cout << " For " << desc[i] << " numbers (~" << scientific << value
         << "): spacing =" << spacing << endl;
  }
}
void machine_epsilon() {
  cout << "\n\033[1m\033[37m\033[4m2. MACHINE EPSILON\033[0m\n\n";
  double c_values[] = {1.0, 10.0, 100.0, 1000.0, 10000.0};
  int num_c = 5;
  cout << "\033[4m\033[1mC" << setw(10) << " | " << setw(10)
       << "Epsilon (theor)"
       << " | "
       << "Epsilon (float)" << " | " << setw(20) << "Epsilon (double)" << " | "
       << setw(20) << "Iterations\033[0m" << endl;
  double theoretical_eps = pow(2.0, -T);
  for (int i = 0; i < num_c; i++) {
    double c = c_values[i];

    float eps_f = 1.0f;
    int iter_f = 0;
    while (static_cast<float>(c) + eps_f > static_cast<float>(c) &&
           iter_f < 1000) {
      eps_f /= 2.0f;
      iter_f++;
    }
    eps_f *= 2.0f;

    double eps_d = 1.0;
    int iter_d = 0;
    while (c + eps_d > c && iter_d < 1000) {
      eps_d /= 2.0;
      iter_d++;
    }
    eps_d *= 2.0;

    cout << scientific << setprecision(2);

    cout << c << " | " << setw(15) << theoretical_eps << " | " << setw(15)
         << eps_f << " | " << setw(20) << eps_d << " | " << setw(15) << iter_d
         << endl;
  }
}
void sum_errors() {
  cout << "\n\033[1m\033[4m\033[37m3. SUM ERRORS\033[0m" << endl << endl;
  double h_values[] = {0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001};
  int num_h = 6;
  cout << "\033[1m\033[4m" << setw(12) << "Step h" << " | " << setw(15)
       << "Sum (float)"
       << " | " << setw(15) << "Sum (double)" << " | " << setw(15)
       << "Exact approx"
       << " | " << setw(12) << "Abs Error" << " | " << setw(12)
       << "Rel Error\033[0m" << endl;
  double ref_sum = 0.0;
  double refh = 1e-6;
  int ref_steps = static_cast<int>(1.0 / refh);
  for (int i = 1; i <= ref_steps; i++) {
    double x = 1.0 + i * refh;
    if (x <= 2.0)
      ref_sum += f(x) * refh;
  }

  for (int i = 0; i < num_h; i++) {
    double h = h_values[i];
    int steps = static_cast<int>(1.0 / h);

    // float sum
    float sum_f = 0.0f;
    for (int j = 1; j <= steps; j++) {
      double x = 1.0f + j * (h);
      if (x <= 2.0f) {
        float f_val = static_cast<float>(f(x));
        sum_f += f_val * static_cast<float>(h);
      }
    }
    // double sum
    double sum_d = 0.0;
    for (int j = 1; j <= steps; j++) {
      double x = 1.0 + j * h;
      if (x <= 2.0) {
        sum_d += f(x) * h;
      }
    }

    double abs_error = fabs(sum_d - ref_sum);
    double rel_error = abs_error / fabs(ref_sum);

    cout << scientific << setprecision(2);
    cout << setw(12) << h << " | " << setw(15) << sum_f << " | " << setw(15)
         << sum_d << " | " << setw(15) << ref_sum << " | " << setw(12)
         << abs_error << " | " << setw(12) << rel_error << endl;
  }
}
double exp_func_dir(double x, double eps, int &terms) {
  double sum = 0.0;
  double term = 1.0;
  terms = 0;

  while (fabs(term) > eps && terms < 100) {
    sum += term;
    terms++;
    term *= x / terms;
  }
  return sum;
}
double exp_func_inv(double x, double eps, int &terms) {
  if (x >= 0) {
    double sum_neg = 0.0;
    double term = 1.0;
    terms = 0;
    while (fabs(term) > eps && terms < 100) {
      sum_neg += term;
      terms++;
      term *= (-x) / terms;
    }
    return 1.0 / sum_neg;
  } else {
    return exp_func_dir(x, eps, terms);
  }
}
void exponential() {
  cout << "\n\033[4m\033[1m4. EXPONENTIAL FUNC\033[0m" << endl << endl;
  double x_vals[] = {-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0};
  int num_x = 7;

  cout << "Comparison of 2 methods of calculating exp(x) with eps=" << EPSILON
       << endl;
  cout << "\033[4m\033[1m" << setw(10) << "x" << " | " << setw(15) << "Exact"
       << " | " << setw(15) << "Method 1" << " | " << setw(10) << "Terms 1"
       << " | " << setw(15) << "Method 2" << " | " << setw(10) << "Terms 2"
       << " | " << setw(12) << "Error 1" << " | " << setw(12) << "Error 2"
       << "\033[0m" << endl;
  for (int i = 0; i < num_x; i++) {
    double x = x_vals[i];
    double exact = exp(x);

    int terms1, terms2;
    double m1 = exp_func_dir(x, EPSILON, terms1);
    double m2 = exp_func_inv(x, EPSILON, terms2);

    double error1 = fabs(m1 - exact);
    double error2 = fabs(m2 - exact);

    cout << scientific << setprecision(2);
    cout << setw(10) << x << " | " << setw(15) << exact << " | " << setw(15)
         << m1 << " | " << setw(10) << terms1 << " | " << setw(15) << m2
         << " | " << setw(10) << terms2 << " | " << setw(12) << error1 << " | "
         << setw(12) << error2 << endl;
  }
}
int main(void) {

  distribution();    // 1
  machine_epsilon(); // 2
  sum_errors();      // 3
  exponential();     // 4;
  return 0;
}
