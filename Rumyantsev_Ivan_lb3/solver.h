#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <string>
using namespace std;

double norm(const vector<double>& v);
vector<double> residual(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b);
bool jacobiMethod(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x, int& iterations,
                  double eps = 1e-8, int maxIter = 50000);
void runTestJacobi(const vector<vector<double>>& A, const vector<double>& b, const string& label);

#endif