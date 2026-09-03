#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <string>
using namespace std;

double norm(const vector<double>& v);
vector<double> residual(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b);
bool gauss(vector<vector<double>> A, vector<double> b, vector<double>& x);
void runTest(const vector<vector<double>>& A, const vector<double>& b, const string& label);

#endif