import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Конфигурация методов
methods = {
    1: {"name": "Bisection Method", "color": "blue", "marker": "o"},
    2: {"name": "Chord Method", "color": "red", "marker": "s"},
    3: {"name": "Newton's Method", "color": "green", "marker": "^"},
    4: {"name": "Simple Iteration Method", "color": "orange", "marker": "d"}
}

def load_data(method_id):
    fname = f"method{method_id}_research.csv"
    if not os.path.exists(fname):
        print(f"Warning: {fname} not found!")
        return None
    # Читаем CSV, который генерирует research.cpp
    return pd.read_csv(fname, sep=';')

def plot_all():
    all_data = {}
    for mid in methods:
        data = load_data(mid)
        if data is not None:
            all_data[mid] = data

    if not all_data:
        print("No data files found. Please run ./research first.")
        return

    # 1. Построение индивидуальных графиков
    for mid, info in methods.items():
        if mid in all_data:
            plt.figure(figsize=(10, 6))
            d = all_data[mid]
            
            plt.semilogx(d['eps'], d['iterations'], marker=info['marker'], 
                         color=info['color'], linewidth=2, markersize=8, label=info['name'])
            
            # Для бисекции теоретическая линия
            if mid == 1:
                a, b = 0.0, 1.0
                eps_theory = np.logspace(np.log10(d['eps'].min()), np.log10(d['eps'].max()), 100)
                n_theory = np.log2((b - a) / eps_theory)
                plt.semilogx(eps_theory, n_theory, '--', color='gray', label='Theoretical N ≈ log₂((b-a)/ε)')

            plt.xlabel('Eps (accuracy)', fontsize=12)
            plt.ylabel('Number of iterations N', fontsize=12)
            plt.title(f"{info['name']}: iterations vs accuracy", fontsize=14)
            plt.grid(True, which='both', linestyle='--', alpha=0.7)
            plt.legend(fontsize=11)
            
            filename = f"{info['name'].lower().replace(' ', '_')}_convergence.png"
            plt.tight_layout()
            plt.savefig(filename, dpi=150)
            print(f"Saved: {filename}")
            plt.close()

    # 2. Общий сравнительный график
    plt.figure(figsize=(10, 6))
    for mid, info in methods.items():
        if mid in all_data:
            d = all_data[mid]
            plt.semilogx(d['eps'], d['iterations'], marker=info['marker'], 
                         color=info['color'], linewidth=2, markersize=8, label=info['name'])

    plt.xlabel('Eps (accuracy)', fontsize=12)
    plt.ylabel('Number of iterations N', fontsize=12)
    plt.title('Comparison of all methods: iterations vs accuracy', fontsize=14)
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.legend(fontsize=11)
    
    plt.tight_layout()
    plt.savefig('comparison_convergence.png', dpi=150)
    print("Saved: comparison_convergence.png")
    plt.show()

if __name__ == "__main__":
    plot_all()
