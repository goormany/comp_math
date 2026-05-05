#include <iostream>
#include <vector>
#include <string>
#include "matrix.h"
#include "solver.h"

using namespace std;

int main() {
    cout << "========================================\n";
    cout << "  Lab: Jacobi Method for SLAE\n";
    cout << "  Student: Rumyantsev I.V., gr. 4343\n";
    cout << "========================================\n\n";

    cout << "Mode:\n";
    cout << "  1 - Manual (choose size)\n";
    cout << "  2 - Auto (tests: 100, 500, 1000)\n";
    cout << "Your choice: ";
    int mode;
    cin >> mode;

    if (mode == 1) {
        cout << "Enter matrix size n: ";
        int n;
        cin >> n;

        cout << "Generating... " << flush;
        auto A = diagonallyDominantMatrix(n);
        auto b = generateB(A);
        cout << "done\n";
        cout << "Diag. dominance: " << (checkDiagonalDominance(A) ? "YES" : "NO") << "\n\n";

        runTestJacobi(A, b, "Manual (" + to_string(n) + "x" + to_string(n) + ")");
    }
    else if (mode == 2) {
        int sizes[] = {100, 500, 1000};
        string names[] = {"DiagDominant 100", "DiagDominant 500", "DiagDominant 1000"};

        for (int t = 0; t < 3; ++t) {
            cout << "[" << t+1 << "/3] " << names[t] << "\n";
            cout << "  Generating... " << flush;
            auto A = diagonallyDominantMatrix(sizes[t]);
            auto b = generateB(A);
            cout << "done (diag.dom: " << (checkDiagonalDominance(A) ? "YES" : "NO") << ")\n";

            runTestJacobi(A, b, names[t]);
        }
    }
    else {
        cerr << "Invalid mode\n";
        return 1;
    }
    return 0;
}