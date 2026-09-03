import csv
import matplotlib.pyplot as plt
import os

methods = {
    "bisection": "Bisection",
    "chord": "Chord",
    "newton": "Newton",
    "iteration": "Simple iteration"
}

for key, label in methods.items():
    fname = f"{key}_iterations_vs_eps.csv"
    if not os.path.exists(fname):
        print(f"File not found: {fname}")
        continue

    eps_values = []
    iterations_values = []

    with open(fname, newline='') as file:
        reader = csv.reader(file, delimiter=';')
        next(reader)
        for row in reader:
            if len(row) < 2:
                continue
            try:
                eps = float(row[0])
                iterations = int(float(row[1]))
            except:
                continue
            eps_values.append(eps)
            iterations_values.append(iterations)

    if len(eps_values) == 0:
        print("No data in", fname)
        continue

    pairs = sorted(zip(eps_values, iterations_values), key=lambda p: p[0])
    eps_sorted = [p[0] for p in pairs]
    it_sorted = [p[1] for p in pairs]

    plt.figure(figsize=(6, 4))
    plt.plot(eps_sorted, it_sorted, marker='o', linestyle='-')
    plt.xscale('log')
    plt.xlabel("Eps")
    plt.ylabel("Iterations")
    plt.title(f"{label}: iterations vs Eps")
    plt.grid(True, which="both", ls="--", alpha=0.6)
    out_png = f"{key}_N_vs_eps.png"
    plt.savefig(out_png, dpi=200, bbox_inches='tight')
    print("Saved plot:", out_png)
    plt.show()

# Вывод CSV в консоль для отчёта
print("\n" + "="*50)
print("CSV DATA FOR REPORT")
print("="*50)

for key, label in methods.items():
    fname_eps = f"{key}_iterations_vs_eps.csv"
    fname_delta = f"{key}_sensitivity_delta.csv"

    print(f"\n--- {label} ---")

    if os.path.exists(fname_eps):
        print(f"\nN(eps):")
        with open(fname_eps, 'r') as f:
            print(f.read())

    if os.path.exists(fname_delta):
        print(f"\nDelta sensitivity (eps=1e-6):")
        with open(fname_delta, 'r') as f:
            print(f.read())