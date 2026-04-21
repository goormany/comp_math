#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

double randDouble(double min = -5.0, double max = 5.0) {
  return min + (max - min) * (rand() / (double)RAND_MAX);
}

void generateDD(int n, ofstream &file) {
  vector<vector<double>> A(n, vector<double>(n));
  vector<double> b(n);

  for (int i = 0; i < n; i++) {
    double rowSum = 0.0;

    for (int j = 0; j < n; j++) {
      if (i != j) {
        A[i][j] = randDouble();
        rowSum += abs(A[i][j]);
      }
    }

    A[i][i] = rowSum + randDouble(1.0, 5.0);
    b[i] = randDouble();
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      file << A[i][j] << " ";
    file << b[i] << "\n";
  }
}

void generateNotDD(int n, ofstream &file) {
  vector<vector<double>> A(n, vector<double>(n));
  vector<double> b(n);

  for (int i = 0; i < n; i++) {
    double rowSum = 0.0;

    for (int j = 0; j < n; j++) {
      if (i != j) {
        A[i][j] = randDouble();
        rowSum += abs(A[i][j]);
      }
    }

    A[i][i] = rowSum * 0.5;
    b[i] = randDouble();
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      file << A[i][j] << " ";
    file << b[i] << "\n";
  }
}

void generateNoSolution(int n, ofstream &file) {
  vector<vector<double>> A(n, vector<double>(n, 0.0));
  vector<double> b(n);

  // Two identical rows with different RHS
  for (int j = 0; j < n; j++) {
    double val = randDouble();
    A[0][j] = val;
    A[1][j] = val;
  }

  b[0] = randDouble();
  b[1] = b[0] + 1.0;

  for (int i = 2; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A[i][j] = randDouble();
    }
    b[i] = randDouble();
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      file << A[i][j] << " ";
    file << b[i] << "\n";
  }
}
void generateSolvable(int n, ofstream &file) {
  vector<vector<double>> A(n, vector<double>(n));
  vector<double> x_true(n);
  vector<double> b(n, 0.0);

  // 1. random matrix A
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A[i][j] = randDouble();
    }
  }

  cout << "\nSolutions:\n";
  for (int i = 0; i < n; i++) {
    x_true[i] = randDouble();
    cout << "x_" << i + 1 << " = " << x_true[i] << "\n";
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      b[i] += A[i][j] * x_true[j];
    }
  }

  // 4. write to file
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      file << A[i][j] << " ";
    }
    file << b[i] << "\n";
  }
}
void generateSolvableDD(int n, ofstream &file) {
  vector<vector<double>> A(n, vector<double>(n));
  vector<double> x_true(n);
  vector<double> b(n, 0.0);

  // 1. choose true solution
  cout << "\nTrue solution (x_true):\n";
  for (int i = 0; i < n; i++) {
    x_true[i] = randDouble();
    cout << "x_" << i + 1 << " = " << x_true[i] << "\n";
  }

  // 2. build diagonally dominant matrix A
  for (int i = 0; i < n; i++) {
    double rowSum = 0.0;

    for (int j = 0; j < n; j++) {
      if (i != j) {
        A[i][j] = randDouble();
        rowSum += fabs(A[i][j]);
      }
    }

    // enforce diagonal dominance
    A[i][i] = rowSum + randDouble(1.0, 5.0);
  }

  // 3. compute b = A * x_true
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      b[i] += A[i][j] * x_true[j];
    }
  }

  // 4. write system to file
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      file << A[i][j] << " ";
    }
    file << b[i] << "\n";
  }
}
int main(int argc, char **argv) {
  srand(time(0));

  string outpath = "";
  string type = "";
  int size = -1;

  for (int i = 1; i < argc; i++) {
    string arg = argv[i];

    if (arg == "--outpath" && i + 1 < argc) {
      outpath = argv[++i];
    } else if (arg == "--size" && i + 1 < argc) {
      size = stoi(argv[++i]);
    } else if (arg == "--type" && i + 1 < argc) {
      type = argv[++i];
    }
  }

  if (outpath.empty() || size <= 0 || type.empty()) {
    cout << "Usage:\n"
         << "./matrix_generator --outpath file.txt --size N --type "
            "dd|not_dd|no_solution\n";
    return 1;
  }

  ofstream file(outpath);
  if (!file) {
    cerr << "Error: cannot open output file\n";
    return 1;
  }

  if (type == "dd") {
    generateDD(size, file);
  } else if (type == "not_dd") {
    generateNotDD(size, file);
  } else if (type == "no_solution") {
    generateNoSolution(size, file);
  } else if(type == "solvable"){
      generateSolvableDD(size,file);
  } else {
    cerr << "Unknown type: " << type << "\n";
    return 1;
  }

  cout << "Generated " << type << " matrix (" << size << "x" << size << ") -> "
       << outpath << "\n";

  return 0;
}
