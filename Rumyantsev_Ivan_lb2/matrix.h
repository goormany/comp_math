#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
using namespace std;

vector<vector<double>> hilbertMatrix(int n);
vector<vector<double>> randomMatrix(int n);
vector<vector<double>> tridiagonalMatrix(int n);
vector<double> generateB(const vector<vector<double>>& A);

#endif