#include <cctype>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct PreciseTimer {
  std::chrono::time_point<std::chrono::high_resolution_clock> start;
  double elapsed = 0.0;

  PreciseTimer() { start = std::chrono::high_resolution_clock::now(); }

  void pause() {
    auto now = std::chrono::high_resolution_clock::now();
    elapsed += std::chrono::duration<double>(now - start).count();
  }

  void resume() { start = std::chrono::high_resolution_clock::now(); }
};

enum class SolutionType { UNIQUE, INFINITE, NONE };

void printMatrix(const std::vector<std::vector<double>> &matrix) {
  if (matrix.empty())
    return;
  int n = matrix.size();
  int cols = matrix[0].size();
  for (int i = 0; i < n; ++i) {
    std::cout << "  [ ";
    for (int j = 0; j < cols; ++j) {
      std::cout << std::setw(9) << std::fixed << std::setprecision(4)
                << matrix[i][j];
      if (j == cols - 2)
        std::cout << " |";
      else if (j < cols - 1)
        std::cout << "  ";
    }
    std::cout << " ]\n";
  }
  std::cout << "\n";
}

SolutionType solveGaussian(std::vector<std::vector<double>> &matrix,
                           std::vector<double> &solution,
                           double &computation_time) {
  computation_time = 0.0;
  PreciseTimer timer;

  int n = matrix.size();
  if (n == 0 || matrix[0].size() != static_cast<size_t>(n + 1)) {
    std::cerr << "Error: Matrix must be n x (n+1).\n";
    return SolutionType::NONE;
  }

  constexpr double TOL = 1e-12;
  int pivot_row = 0;
  int step = 1;

  timer.pause();
  std::cout << "\033[33m[Step " << step++
            << "]\033[0m Initial augmented matrix:\n";
  printMatrix(matrix);
  timer.resume();

  // transform to row-echelon form using forward elemination
  for (int col = 0; col < n && pivot_row < n; ++col) {
    // partially pivot to find row with largest |value| in current column
    int max_row = pivot_row;
    double max_val = std::abs(matrix[pivot_row][col]);
    for (int i = pivot_row + 1; i < n; ++i) {
      if (std::abs(matrix[i][col]) > max_val) {
        max_val = std::abs(matrix[i][col]);
        max_row = i;
      }
    }

    // If entire column below pivot_row is ~0, this variable is dependent/free
    if (max_val < TOL) {
      timer.pause();
      std::cout << "\033[33m[Step " << step++ << "]\033[0m Column " << (col + 1)
                << " has no significant pivot. Skipping (indicates free "
                   "variable).\n\n";
      timer.resume();
      continue;
    }

    // Swap the pivot row into position
    if (max_row != pivot_row) {
      timer.pause();
      std::cout << "\033[33m[Step " << step++ << "]\033[0m Swap Row "
                << (pivot_row + 1) << " and Row " << (max_row + 1) << "\n";
      timer.resume();

      std::swap(matrix[pivot_row], matrix[max_row]);

      timer.pause();
      printMatrix(matrix);
      timer.resume();
    }

    // Eliminate entries below the pivot
    for (int i = pivot_row + 1; i < n; ++i) {
      double factor = matrix[i][col] / matrix[pivot_row][col];

      timer.pause();
      std::cout << "\033[33m[Step " << step++ << "]\033[0m Row " << (i + 1)
                << " = Row " << (i + 1) << " - (" << factor << ") * Row "
                << (pivot_row + 1) << "\n";
      timer.resume();

      for (int j = col; j <= n; ++j) {
        matrix[i][j] -= factor * matrix[pivot_row][j];
      }

      timer.pause();
      printMatrix(matrix);
      timer.resume();
    }
    pivot_row++;
  }

  int rank = pivot_row; // Number of independent equations

  // 0x + 0y + ... = non-zero
  for (int i = rank; i < n; ++i) {
    if (std::abs(matrix[i][n]) > TOL) {
      timer.pause();
      std::cout << "\033[33m[Step " << step++
                << "]\033[0m Inconsistent matrix: Row " << (i + 1)
                << " implies 0 = " << matrix[i][n] << "\n\n";
      timer.resume();
      return SolutionType::NONE;
    }
  }

  // Check if has infinite
  if (rank < n) {
    timer.pause();
    std::cout << "\033[33m[Step " << step++
              << "]\033[0m System is consistent but underdetermined. "
              << (n - rank) << " free variable(s) exist.\n\n";
    timer.resume();
    return SolutionType::INFINITE;
  }

  // Back substitution
  // x_i = (b_i-S(a_ij * x_i))/a_ii
  // Take the rhs b_i and subtract all the already known
  // vars from x[i+1] to x[n-1] and divide by the pivot
  timer.pause();
  std::cout << "\033[33m[Step " << step++
            << "]\033[0m Matrix is in upper triangular form.\n";
  printMatrix(matrix);
  timer.resume();

  solution.resize(n);
  for (int i = n - 1; i >= 0; --i) {
    double rhs = matrix[i][n];
    double subtracted = 0.0;

    timer.pause();
    std::cout << "\033[33m[Step " << step++
              << "]\033[0m Back-substitution: Calculate x" << (i + 1) << "\n";
    std::cout << "  x" << (i + 1) << " = ( " << rhs;

    for (int j = i + 1; j < n; ++j) {
      double term = matrix[i][j] * solution[j];
      subtracted += term;
      std::cout << " - (" << matrix[i][j] << " * " << solution[j] << ")";
    }

    std::cout << " ) / " << matrix[i][i];
    timer.resume();

    // Actual computation happens while timer is running
    solution[i] = (rhs - subtracted) / matrix[i][i];

    timer.pause();
    std::cout << " = " << solution[i] << "\n\n";
    timer.resume();
  }

  computation_time = timer.elapsed;
  return SolutionType::UNIQUE;
}

const char *solutionTypeToString(SolutionType type) {
  switch (type) {
  case SolutionType::UNIQUE:
    return "UNIQUE";
  case SolutionType::INFINITE:
    return "INFINITE";
  case SolutionType::NONE:
    return "NO SOLUTION";
  default:
    return "UNKNOWN";
  }
}

bool readMatrixFromFile(const std::string &filename,
                        std::vector<std::vector<double>> &matrix) {
  std::ifstream in(filename);
  if (!in.is_open()) {
    std::cerr << "Error: Could not open file '" << filename << "'.\n";
    return false;
  }

  matrix.clear();
  std::vector<double> flat;
  std::string line;

  while (std::getline(in, line)) {
    bool is_blank = true;
    for (char c : line) {
      if (!std::isspace(static_cast<unsigned char>(c))) {
        is_blank = false;
        break;
      }
    }
    if (is_blank)
      continue;

    std::istringstream iss(line);
    double val;
    while (iss >> val) {
      flat.push_back(val);
    }
  }
  in.close();

  int total = flat.size();
  if (total == 0) {
    std::cerr << "Error: File is empty or contains no numbers.\n";
    return false;
  }

  double discriminant = std::sqrt(1.0 + 4.0 * total);
  int n = static_cast<int>((-1.0 + discriminant) / 2.0);

  if (n <= 0 || n * (n + 1) != total) {
    std::cerr << "Error: File contains " << total
              << " values, which cannot form an n x (n+1) augmented matrix.\n";
    return false;
  }

  matrix.resize(n);
  int idx = 0;
  for (int i = 0; i < n; ++i) {
    matrix[i].reserve(n + 1);
    for (int j = 0; j <= n; ++j) {
      matrix[i].push_back(flat[idx++]);
    }
  }
  return true;
}

int main(int argc, char *argv[]) {

  if (argc != 2) {
    std::cerr << "\033[31mUsage:\033[0m " << argv[0] << " <matrix_file.txt>\n";
    return 1;
  }

  std::string filename = argv[1];
  std::vector<std::vector<double>> matrix;

  if (!readMatrixFromFile(filename, matrix)) {
    return 1;
  }

  std::vector<double> solution;
  double compute_time = 0.0;
  auto *orig = std::cout.rdbuf();
  std::cout.rdbuf(NULL);
  SolutionType result = solveGaussian(matrix, solution, compute_time);
  std::cout.rdbuf(orig);

  std::cout << "\033[4m\033[1mFinal Result\033[0m\n";
  std::cout << "Solution type: " << solutionTypeToString(result) << "\n";
  std::cout << "Computation time (excluding I/O): " << std::fixed
            << std::setprecision(6) << compute_time << " seconds\n";

  std::cout << std::fixed << std::setprecision(10);

  if (result == SolutionType::UNIQUE) {
    std::cout << "\nSolution vector x:\n";
    for (int i = 0; i < static_cast<int>(matrix.size()); ++i) {
      std::cout << "  x[" << i << "] = " << solution[i] << "\n";
    }
  } else if (result == SolutionType::INFINITE) {
    std::cout
        << "The system has infinitely many solutions (free variables exist).\n";
  } else {
    std::cout << "The system is inconsistent (no solution exists).\n";
  }

  return 0;
}
