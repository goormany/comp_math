#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#define EPS 1e-9
using namespace std;

bool isDD(const vector<vector<double>> &A) {
  int n = A.size();

  for (int i = 0; i < n; i++) {
    double diag = fabs(A[i][i]);
    double rowSum = 0.0;

    for (int j = 0; j < n; j++) {
      if (j != i) {
        rowSum += fabs(A[i][j]);
      }
    }
    if (diag < rowSum)
      return false;
  }
  return true;
}
double computeResidual(const vector<vector<double>> &A,
                       const vector<double> &x,
                       const vector<double> &b) {
  int n = A.size();
  double err = 0.0;

  for (int i = 0; i < n; i++) {
    double Ax_i = 0.0;

    for (int j = 0; j < n; j++) {
      Ax_i += A[i][j] * x[j];
    }

    err += fabs(Ax_i - b[i]);
  }

  return err;
}

vector<double> jacobiSolver(const vector<vector<double>> &A,
                            const vector<double> &b, int maxIter = 5000,
                            double eps = EPS) {
  int n = A.size();
  vector<double> x(n, 0.0), x_new(n, 0.0);
  double min = 1e-12;

  for (int iter = 0; iter < maxIter; iter++) {
    for (int i = 0; i < n; i++) {
      double sum = 0.0;

      for (int j = 0; j < n; j++) {
        if (j != i) {
          sum += A[i][j] * x[j];
        }
      }

      if (fabs(A[i][i]) < min) {
        throw runtime_error("0 diagonal element");
      }

      x_new[i] = (b[i] - sum) / A[i][i];
    }

    double maxErr = 0.0;
    for (int i = 0; i < n; i++) {
      maxErr = max(maxErr, fabs(x_new[i] - x[i]));
    }

    if (maxErr < eps) {
      cout << "Converged in " << iter + 1 << " Iterations.\n";
      return x_new;
    }
    x = x_new;
  }
  cout << "Max iterations reached without convergence." << endl;
  return x;
}

void readFromFile(const string &filename, vector<vector<double>> &A,
                  vector<double> &b) {

  ifstream file(filename);

  if (!file) {
    throw runtime_error("Cannot open file");
  }

  vector<vector<double>> tmp;
  double value;

  while (true) {
    vector<double> row;
    for (int i = 0;; i++) {
      if (!(file >> value))
        break;
      row.push_back(value);

      if (file.peek() == '\n')
        break;
    }
    if (row.empty())
      break;
    tmp.push_back(row);
  }

  int n = tmp.size();

  for (const auto &row : tmp) {
    if ((int)row.size() != n + 1) {
      throw runtime_error("Invalid matrix size.");
    }
  }

  A.resize(n, vector<double>(n));
  b.resize(n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A[i][j] = tmp[i][j];
    }
    b[i] = tmp[i][n];
  }
}

int main(int argc, char **argv) {
  if (argc != 2) {
    cout << "\033[31mError: \033[0m" << "Usage: ./[main] file.txt" << std::endl;
    exit(1);
  }
  string filename = argv[1];
  vector<vector<double>> A;
  vector<double> b;

  try {
    readFromFile(filename, A, b);
    int n = A.size();
    cout << "Matrix of size " << n << "x" << n << " loaded" << endl;
    if (!isDD(A)) {
      cout << "Matrix not diagonally dominant" << endl;
      return 1;
    } else
      cout << "Diagonally dominant matrix" << endl;

    vector<double> sln = jacobiSolver(A, b);
    double residual = computeResidual(A, sln, b);
    double b_norm = 0.0;
    for (double v : b) b_norm += fabs(v);
    double rres = (b_norm == 0) ? residual : residual/b_norm;

    cout << "\n\033[1m\033[4mSolution check: \033[0m\n";
    cout << "||Ax-b|| = " << residual << endl;
    cout << "Relative residual = " << rres << endl;
    if(rres < EPS){
        cout << "Good solution" << endl;
    }else if(rres < 1e-6){
        cout << "Approx solution" << endl;
    }else{
        cout << "Bad solution" << endl;
    }
    cout << "\033[1m\033[4mSolutions:\033[0m" << endl;
    int num = 0;
    for (double x : sln) {
      num++;
      //cout << "x_" << num << " = " << x << endl;
    }
    puts("");
  } catch (const exception &e) {
    cerr << "\033[31mError: " << e.what() << endl;
  }
  return 1;
}
