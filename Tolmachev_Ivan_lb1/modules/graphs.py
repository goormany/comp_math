import csv
import matplotlib.pyplot as plt
import os
import math

methods = {
    1: "Бисекция",
    2: "Хорд",
    3: "Ньютон",
    4: "Простые итерации"
}

for method in range(1,5):
    fname = f"method{method}_iterations_vs_eps.csv"
    if not os.path.exists(fname):
        print(f"Файл не найден: {fname}. Сначала запустите C++ программу.")
        continue

    eps_values = []
    iterations_values = []
    roots = []

    with open(fname, newline='') as file:
        reader = csv.reader(file, delimiter=';')
        header = next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            try:
                eps = float(row[0])
                iterations = int(float(row[1]))
                root = float(row[2]) if len(row) > 2 else float('nan')
            except:
                continue
            eps_values.append(eps)
            iterations_values.append(iterations)
            roots.append(root)

    if len(eps_values) == 0:
        print("Нет данных в", fname)
        continue

    pairs = sorted(zip(eps_values, iterations_values), key=lambda p: p[0])
    eps_sorted = [p[0] for p in pairs]
    it_sorted = [p[1] for p in pairs]

    plt.figure(figsize=(6,4))
    plt.plot(eps_sorted, it_sorted, marker='o', linestyle='-')
    plt.xscale('log')
    plt.xlabel("Eps")
    plt.ylabel("Iterations")
    plt.title(f"{methods[method]}: число итераций vs Eps")
    plt.grid(True, which="both", ls="--", alpha=0.6)
    out_png = f"method{method}_N_vs_eps.png"
    plt.savefig(out_png, dpi=200, bbox_inches='tight')
    print("Saved plot:", out_png)
    plt.show()