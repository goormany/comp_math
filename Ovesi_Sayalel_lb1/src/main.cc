#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

// f(x)=2^(x^2)-1/x
double F(double x) { return pow(2.0, x * x) - 1.0 / x; }

// derivative: f'(x)=2^(x^2)*ln(2)*2x+1/(x^2)
double FDeriv(double x) {
  double pow_term = pow(2.0, x * x);
  return pow_term * log(2.0) * 2.0 * x + 1.0 / (x * x);
}

double Round(double value, int decs) {
  if (decs < 0)
    return value;
  double fact = pow(10.0, decs);
  return floor(value * fact) / fact;
}

double BISECT(double Left, double Right, double Eps, int &N) {
  double FLeft = F(Left);
  double FRight = F(Right);
  double x;
  if (FLeft * FRight > 0.0) {
    cout << "\033[31mError:\033[0m Function has same sign at interval ends"
         << endl;
    exit(1);
  }
  if (Eps <= 0.0) {
    cout << "\033[31mError:\033[0m Invalid precision" << endl;
    exit(1);
  }
  N = 0;
  while ((Right - Left) > 2 * Eps) {
    x = (Left + Right) / 2.0;
    double Fx = F(x);

    if (Fx == 0.0)
      return x;
    if (FLeft * Fx < 0.0) {
      Right = x;
      FRight = Fx;
    } else {
      Left = x;
      FLeft = Fx;
    }
    N++;
  }
  return (Left + Right) / 2.0;
}

double HORDA(double Left, double Right, double Eps, int &N) {
  double FLeft = F(Left);
  double FRight = F(Right);
  double X, Y;
  if (FLeft * FRight > 0.0) {
    cout << "\033[31mError:\033[0m Invalid interval" << endl;
    exit(1);
  }
  if (Eps <= 0.0) {
    cout << "\033[31mError:\033[0m Invalid precision" << endl;
    exit(1);
  }
  if (FLeft == 0.0)
    return Left;
  if (FRight == 0.0)
    return Right;

  N = 0;
  do {
    X = Left - (Right - Left) * FLeft / (FRight - FLeft);
    Y = F(X);
    if (Y == 0.0)
      return X;

    if (Y * FLeft < 0.0) {
      Right = X;
      FRight = Y;
    } else {
      Left = X;
      FLeft = Y;
    }
    N++;
  } while (fabs(Y) >= Eps);

  return X;
}

void stabanas_bisect(double Left, double Right, int maxIter) {
  cout << "\n\033[4m\033[1m\033[37mBISECT STABILITY ANALYSIS\033[0m" << endl;

  // Find accurate root
  int N_acc;
  double x_acc = BISECT(Left, Right, 1e-15, N_acc);
  cout << "Accurate root: " << setprecision(15) << x_acc << endl;
  cout << "Iterations for accurate root: " << N_acc << endl;

  double nu_delta = 1.0 / fabs(FDeriv(x_acc));
  cout << "nu_delta (condition number): " << nu_delta << endl;

  vector<int> decimals = {6, 5, 4, 3, 2, 1};
  cout << "\nAnalysis for different perturbation levels:\n";
  vector<double> first_stable_iter;

  for (int dec : decimals) {
    double delta = pow(10.0, -dec);
    cout << "\n--- Perturbation to " << dec
         << " decimal places (delta = " << delta << ") ---\n";

    vector<double> delta_i; // Differences |x_N - x_i*|

    double left = Left, right = Right;
    double FLeft = Round(F(left), dec);
    double FRight = Round(F(right), dec);

    int iter = 0;
    int first_stable = -1;

    cout << "+-----+---------------+---------------+---------------+-----------"
            "----+---------------+\n";
    cout << "|" << setw(4) << "i" << " |" << setw(14) << "x_i" << " |"
         << setw(14) << "x_i*"
         << " |" << setw(14) << "Δ_i" << " |" << setw(14) << "ν_Δ*ε" << " |" << 
         setw(14) << "Condition" << " |\n";
    cout << "+-----+---------------+---------------+---------------+-----------"
            "----+---------------+\n";

    while (iter < maxIter) {
      // BISECTION step (not chord method!)
      double x = (left + right) / 2.0;

      double Fx_perturbed = Round(F(x), dec);

      // Calculate difference from accurate root
      double diff = fabs(x_acc - x);
      delta_i.push_back(diff);

      double nu_eps = nu_delta * delta;
      string condition =
          (diff <= nu_eps) ? "Well-conditioned" : "Ill-conditioned";

      cout << "|" << setw(4) << iter + 1 << " |" << setw(14) << fixed
           << setprecision(8) << x << " |" << setw(14) << x << " |" << setw(14)
           << scientific << setprecision(4) << diff << " |" << setw(14) << fixed
           << setprecision(8) << nu_eps << " |" << setw(14) << condition
           << " |\n";

      if (diff <= nu_eps && first_stable == -1) {
        first_stable = iter + 1;
      }

      if (Fx_perturbed == 0.0)
        break;

      // Update interval based on perturbed function values
      if (FLeft * Fx_perturbed < 0.0) {
        right = x;
        FRight = Fx_perturbed;
      } else {
        left = x;
        FLeft = Fx_perturbed;
      }

      iter++;
    }
    cout << "+-----+---------------+---------------+---------------+-----------"
            "----+---------------+\n";

    cout << "\nFirst well-conditioned iteration: "
         << (first_stable > 0 ? to_string(first_stable) : "Never stable")
         << endl;

    first_stable_iter.push_back(first_stable);
  }

  // Save data for graph
  ofstream outFile("bisect_stability_data.txt");
  outFile << "Decimals\tDelta\tFirstStableIteration\n";
  for (size_t i = 0; i < decimals.size(); i++) {
    outFile << decimals[i] << "\t" << pow(10.0, -decimals[i]) << "\t"
            << first_stable_iter[i] << "\n";
  }
  outFile.close();
  cout << "\nData saved to bisect_stability_data.txt\n";
}

void stabanas_horda(double Left, double Right, int maxIter) {

  cout << "\n\033[4m\033[1mCHORD METHOD STABILITY ANALYSIS\033[0m\n";

  // First find accurate solution with many iterations
  int N_accurate;
  double x_accurate = HORDA(Left, Right, 1e-15, N_accurate);

  cout << "Accurate root: " << setprecision(15) << x_accurate << endl;
  cout << "Iterations for accurate root: " << N_accurate << endl;

  // Calculate condition number
  double nu_delta = 1.0 / fabs(FDeriv(x_accurate));
  cout << "nu_delta (condition number): " << nu_delta << endl;

  // Analyze for different perturbation levels
  vector<int> decimals = {6, 5, 4, 3, 2, 1};

  cout << "\nAnalysis for different perturbation levels:\n";
  vector<double> first_stable_iter;

  for (int dec : decimals) {
    double delta = pow(10.0, -dec);
    cout << "\n--- Perturbation to " << dec
         << " decimal places (delta = " << delta << ") ---\n";

    vector<double> delta_i; // Differences |x_N - x_i*|

    // Generate sequence with perturbation
    double left = Left, right = Right;
    double FLeft = Round(F(left), dec);
    double FRight = Round(F(right), dec);

    int iter = 0;
    int first_stable = -1;

    cout << "+-----+---------------+---------------+---------------+-----------"
            "----+---------------+\n";
    cout << "|" << setw(4) << "i" << " |" << setw(14) << "x_i" << " |"
         << setw(14) << "x_i*"
         << " |" << setw(14) << "Δ_i" << " |" << setw(14) << "ν_Δ*ε" << " |" << 
         setw(14) << "Condition" << " |\n";
    cout << "+-----+---------------+---------------+---------------+-----------"
            "----+---------------+\n";

    while (iter < maxIter) {
      double x =
          left - (right - left) * FLeft / (FRight - FLeft); // Chord method step

      // Calculate difference from accurate root
      double diff = fabs(x_accurate - x);
      delta_i.push_back(diff);

      double nu_eps = nu_delta * delta;
      string condition =
          (diff <= nu_eps) ? "Well-conditioned" : "Ill-conditioned";

      cout << "|" << setw(4) << iter + 1 << " |" << setw(14) << fixed
           << setprecision(8) << x << " |" << setw(14) << x << " |" << setw(14)
           << scientific << setprecision(4) << diff << " |" << setw(14) << fixed
           << setprecision(8) << nu_eps << " |" << setw(14) << condition
           << " |\n";

      if (diff <= nu_eps && first_stable == -1) {
        first_stable = iter + 1;
      }

      // Compute perturbed value for next iteration
      double Fx_perturbed = Round(F(x), dec);

      if (Fx_perturbed == 0.0)
        break;

      if (Fx_perturbed * FLeft < 0.0) {
        right = x;
        FRight = Fx_perturbed;
      } else {
        left = x;
        FLeft = Fx_perturbed;
      }

      iter++;
    }
    cout << "+-----+---------------+---------------+---------------+-----------"
            "----+---------------+\n";

    cout << "\nFirst well-conditioned iteration: "
         << (first_stable > 0 ? to_string(first_stable) : "Never stable")
         << endl;

    first_stable_iter.push_back(first_stable);
  }

  // Save data for graph
  ofstream outFile("horda_stability_data.txt");
  outFile << "Decimals\tDelta\tFirstStableIteration\n";
  for (size_t i = 0; i < decimals.size(); i++) {
    outFile << decimals[i] << "\t" << pow(10.0, -decimals[i]) << "\t"
            << first_stable_iter[i] << "\n";
  }
  outFile.close();
  cout << "\nData saved to horda_stability_data.txt\n";
}

// Part 1: Study iteration count vs precision for both methods
void IterationVsPrecision(double Left, double Right) {
  vector<double> eps_values = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};

  cout << "\n\033[1m\033[4m\033[37m1. ITERATIONS VS PRECISION\033[0m" << endl;

  cout << "\nBISECTION METHOD:\n";
  cout << "+-------------+---------------+-------------+\n";
  cout << "|" << setw(12) << "Eps" << " |" << setw(14) << "Root" << " |"
       << setw(12) << "Iterations" << " |\n";
  cout << "+-------------+---------------+-------------+\n";
  for (double eps : eps_values) {
    int N;
    double root = BISECT(Left, Right, eps, N);
    cout << "|" << setw(12) << eps << " |" << setw(14) << setprecision(10)
         << root << " |" << setw(12) << N << " |\n";
  }
  cout << "+-------------+---------------+-------------+\n";

  cout << "\nCHORD METHOD (HORDA):\n";
  cout << "+-------------+---------------+-------------+\n";
  cout << "|" << setw(12) << "Eps" << " |" << setw(14) << "Root" << " |"
       << setw(12) << "Iterations" << " |\n";
  cout << "+-------------+---------------+-------------+\n";
  for (double eps : eps_values) {
    int N;
    double root = HORDA(Left, Right, eps, N);
    cout << "|" << setw(12) << eps << " |" << setw(14) << setprecision(10)
         << root << " |" << setw(12) << N << " |\n";
  }
  cout << "+-------------+---------------+-------------+\n";

  // Save data for plotting
  ofstream iterFile("iteration_comparison.txt");
  iterFile << "Eps\tBisection_iters\tHorda_iters\n";
  for (double eps : eps_values) {
    int N_bisect, N_horda;
    BISECT(Left, Right, eps, N_bisect);
    HORDA(Left, Right, eps, N_horda);
    iterFile << eps << "\t" << N_bisect << "\t" << N_horda << "\n";
  }
  iterFile.close();
}

// Function to verify the root interval
bool VerifyInterval(double Left, double Right) {
  cout << "Function f(x) = 2^(x^2) - 1/x\n";
  cout << "Checking interval [" << Left << ", " << Right << "]:\n";
  cout << "f(Left) = " << F(Left) << endl;
  cout << "f(Right) = " << F(Right) << endl;

  if (F(Left) * F(Right) < 0) {
    cout << "-> \033[37mRoot exists in this interval\033[0m\n\n";
    return true;
  } else {
    cout
        << "\033[31mError:\033[0m No root or multiple roots in this interval\n";
    return false;
  }
}

int main() {
  // Find correct interval (avoid x=0 where function is undefined)
  cout << "Searching for root interval...\n";
  cout << "+-------+---------------+\n";
  cout << "|" << setw(6) << "x" << " |" << setw(14) << "f(x)" << " |\n";
  cout << "+-------+---------------+\n";
  for (double x = 0.5; x <= 1.0; x += 0.1) {
    cout << "|" << setw(6) << fixed << setprecision(2) << x << " |" << setw(14)
         << setprecision(6) << F(x) << " |\n";
  }
  cout << "+-------+---------------+\n";

  // Root is between 0.7 and 0.8 (where sign changes)
  double Left = 0.7;
  double Right = 0.8;

  cout << "\nUsing interval [" << Left << ", " << Right << "]\n";
  VerifyInterval(Left, Right);

  // Part 1: Iterations vs Precision
  IterationVsPrecision(Left, Right);

  // Part 2: Stability analysis for both methods
  // cout << "\n" << string(60, '=') << "\n";
  // cout << "PART 2: STABILITY ANALYSIS\n";
  // cout << string(60, '=') << "\n";

  stabanas_bisect(Left, Right, 50);
  stabanas_horda(Left, Right, 30); // Horda converges faster

  // Comparative summary
  // cout << "\n========== COMPARATIVE SUMMARY ==========\n";
  // cout << "1. Both methods successfully found the root of f(x) = 2^(x^2) - "
  // "1/x\n";
  // cout << "2. The chord method (HORDA) converges faster than bisection\n";
  // cout << "3. Stability analysis shows the impact of function value "
  // "perturbations\n";
  // cout << "4. Data files saved for graphing:\n";
  // cout << "   - iteration_comparison.txt\n";
  // cout << "   - bisect_stability_data.txt\n";
  // cout << "   - horda_stability_data.txt\n";

  return 0;
}
