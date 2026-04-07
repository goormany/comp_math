#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>

using namespace std;
float min_value = -17.0;
float max_value = 17.0;
int n;
vector<vector<float>> A;
vector<float> b;


void print_result() {
    int n = A.size();
    cout << n << endl;
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << A[i][j];
            if (j < n - 1) cout << " ";
        }
        cout << endl;
    }
    
    for (int i = 0; i < n; ++i) {
        cout << b[i];
        if (i < n - 1) cout << " ";
    }

    cout << endl;
}


float random_float(float min, float max, mt19937& gen) {
    uniform_real_distribution<float> dist(min, max);
    return dist(gen);
}



void generate(mt19937& gen) {
    A.assign(n, vector<float>(n, 0.0));
    b.assign(n, 0.0);

    uniform_real_distribution<float> gen_element(min_value, max_value);

    for (int i = 0; i < n; i++) {
        float not_diagon_sum = 0.0;

        for (int j = 0; j < n; j++) {
            if (j != i) {
                A[i][j] = gen_element(gen);
                not_diagon_sum += fabs(A[i][j]);
            }
        }

        float diagon_element = not_diagon_sum + random_float(1.0, not_diagon_sum / 3.0, gen);

        if (random_float(0.0, 1.0, gen) < 0.5)
            diagon_element = -diagon_element;

        A[i][i] = diagon_element;
    }

    for (int i = 0; i < n; ++i) {
        b[i] = gen_element(gen);
    }
}



int main() {
    cin >> n;
    
    if (n <= 0) {
        cerr << "Ошибка: n должно быть положительным целым числом.\n";
        return 1;
    }
    
    cin >> min_value;
    cin >> max_value;
    
    if (min_value >= max_value) {
        cerr << "Ошибка: min должно быть меньше max.\n";
        return 1;
    }

    random_device rd;
    mt19937 gen(rd());

 
    
    generate(gen);
    
    print_result();

    return 0;
}
