#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
using namespace std;

vector<vector<double>> diagonallyDominantMatrix(int n);
vector<double> generateB(const vector<vector<double>>& A);
bool checkDiagonalDominance(const vector<vector<double>>& A);

#endif