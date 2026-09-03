#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include "matrix.h"
#include "solver.h"

using namespace std;

void writeMatrix(const vector<vector<double>>& A, const vector<double>& b, const string& filename) {
    int n = A.size();
    ofstream fout(filename, ios::binary);
    fout.write(reinterpret_cast<const char*>(&n), sizeof(int));
    for (int i = 0; i < n; ++i)
        fout.write(reinterpret_cast<const char*>(A[i].data()), n * sizeof(double));
    fout.write(reinterpret_cast<const char*>(b.data()), n * sizeof(double));
    fout.close();
}

void readMatrix(const string& filename, vector<vector<double>>& A, vector<double>& b) {
    ifstream fin(filename, ios::binary);
    int n;
    fin.read(reinterpret_cast<char*>(&n), sizeof(int));
    A.assign(n, vector<double>(n));
    b.resize(n);
    for (int i = 0; i < n; ++i)
        fin.read(reinterpret_cast<char*>(A[i].data()), n * sizeof(double));
    fin.read(reinterpret_cast<char*>(b.data()), n * sizeof(double));
    fin.close();
}

int main() {
    cout << "========================================\n";
    cout << "  Lab #2: Gauss Method for SLAE\n";
    cout << "  Student: Rumyantsev I.V., gr. 4343\n";
    cout << "========================================\n\n";
    cout << "Mode:\n";
    cout << "  1 - Manual (choose type and size)\n";
    cout << "  2 - Auto (run all 4 tests)\n";
    cout << "Your choice: ";

    int mode;
    cin >> mode;

    if (mode == 1) {
        cout << "\nMatrix types:\n";
        cout << "  1 - Hilbert\n";
        cout << "  2 - Random [-1, 1]\n";
        cout << "  3 - Tridiagonal\n";
        cout << "Your choice: ";
        int type;
        cin >> type;

        cout << "Enter size n: ";
        int n;
        cin >> n;

        cout << "Generating... " << flush;
        vector<vector<double>> A;
        switch (type) {
            case 1: A = hilbertMatrix(n); break;
            case 2: A = randomMatrix(n); break;
            case 3: A = tridiagonalMatrix(n); break;
            default: cerr << "Invalid type\n"; return 1;
        }
        vector<double> b = generateB(A);
        cout << "done\n";

        writeMatrix(A, b, "input.bin");

        vector<vector<double>> A_read;
        vector<double> b_read;
        readMatrix("input.bin", A_read, b_read);

        runTest(A_read, b_read, "Manual test");
    }
    else if (mode == 2) {
        string names[] = {"Tridiagonal 1000", "Random 1000", "Hilbert 10", "Hilbert 1000"};
        string files[] = {"tridiagonal_1000.bin", "random_1000.bin",
                         "hilbert_10.bin", "hilbert_1000.bin"};

        cout << "=== Generating ===\n";
        cout << "[1/4] Tridiagonal 1000... " << flush;
        auto A1 = tridiagonalMatrix(1000);
        writeMatrix(A1, generateB(A1), files[0]);
        cout << "done\n";

        cout << "[2/4] Random 1000... " << flush;
        auto A2 = randomMatrix(1000);
        writeMatrix(A2, generateB(A2), files[1]);
        cout << "done\n";

        cout << "[3/4] Hilbert 10... " << flush;
        auto A3 = hilbertMatrix(10);
        writeMatrix(A3, generateB(A3), files[2]);
        cout << "done\n";

        cout << "[4/4] Hilbert 1000... " << flush;
        auto A4 = hilbertMatrix(1000);
        writeMatrix(A4, generateB(A4), files[3]);
        cout << "done\n\n";

        cout << "=== Solving ===\n";
        for (int t = 0; t < 4; ++t) {
            cout << "[" << t+1 << "/4] " << names[t] << "\n";
            vector<vector<double>> A;
            vector<double> b;
            readMatrix(files[t], A, b);
            runTest(A, b, names[t]);
        }
    }
    else {
        cerr << "Invalid mode\n";
        return 1;
    }
    return 0;
}