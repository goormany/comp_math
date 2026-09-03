import csv
import matplotlib.pyplot as plt
import os

methods = [
    ("Бисекция", "bisection", "blue"),
    ("Хорд", "chord", "green"),
    ("Ньютона", "newton", "red"),
    ("Простые итерации", "iteration", "purple")
]

for method_name, prefix, color in methods:
    eps_file = f"{prefix}_eps.csv"
    delta_file = f"{prefix}_delta.csv"

    if not os.path.exists(eps_file):
        print(f"Файл {eps_file} не найден. Пропускаем метод {method_name}.")
    else:
        eps_values = []
        iter_values = []
        roots = []

        with open(eps_file, newline='') as f:
            reader = csv.reader(f, delimiter=';')
            header = next(reader, None)
            for row in reader:
                if len(row) < 3:
                    continue
                try:
                    eps = float(row[0])
                    it = int(row[1])
                    r = float(row[2])
                except ValueError:
                    continue
                eps_values.append(eps)
                iter_values.append(it)
                roots.append(r)

        if eps_values:
            pairs = sorted(zip(eps_values, iter_values), key=lambda p: p[0])
            eps_sorted = [p[0] for p in pairs]
            it_sorted = [p[1] for p in pairs]

            plt.figure(figsize=(8, 5))
            plt.plot(eps_sorted, it_sorted, marker='o', linestyle='-', color=color)
            plt.xscale('log')
            plt.xlabel('Точность ε')
            plt.ylabel('Число итераций')
            plt.title(f'Метод {method_name}: зависимость числа итераций от точности ε')
            plt.grid(True, which='both', linestyle='--', alpha=0.7)
            plt.gca().invert_xaxis()
            plt.tight_layout()
            out_png = f"{prefix}_eps_iter.png"
            plt.savefig(out_png, dpi=200)
            print(f"Сохранён график: {out_png}")
            plt.show()
        else:
            print(f"Нет данных в {eps_file}")

    if not os.path.exists(delta_file):
        print(f"Файл {delta_file} не найден. Пропускаем.")
    else:
        delta_values = []
        iter_values_delta = []
        roots_delta = []

        with open(delta_file, newline='') as f:
            reader = csv.reader(f, delimiter=';')
            header = next(reader, None)
            for row in reader:
                if len(row) < 3:
                    continue
                try:
                    d = float(row[0])
                    it = int(row[1])
                    r = float(row[2])
                except ValueError:
                    continue
                delta_values.append(d)
                iter_values_delta.append(it)
                roots_delta.append(r)

        if delta_values:
            triples = sorted(zip(delta_values, iter_values_delta, roots_delta), key=lambda p: p[0])
            delta_sorted = [t[0] for t in triples]
            it_sorted_delta = [t[1] for t in triples]
            root_sorted_delta = [t[2] for t in triples]

            plt.figure(figsize=(8, 5))
            plt.plot(delta_sorted, root_sorted_delta, marker='s', linestyle='-', color=color)
            plt.xscale('log')
            plt.xlabel('Параметр округления Δ')
            plt.ylabel('Найденный корень')
            plt.title(f'Метод {method_name}: влияние ошибок округления на значение корня')
            plt.grid(True, which='both', linestyle='--', alpha=0.7)
            plt.gca().invert_xaxis()
            plt.tight_layout()
            out_png = f"{prefix}_delta_root.png"
            plt.savefig(out_png, dpi=200)
            print(f"Сохранён график: {out_png}")
            plt.show()

            plt.figure(figsize=(8, 5))
            plt.plot(delta_sorted, it_sorted_delta, marker='^', linestyle='-', color=color)
            plt.xscale('log')
            plt.xlabel('Параметр округления Δ')
            plt.ylabel('Число итераций')
            plt.title(f'Метод {method_name}: влияние ошибок округления на число итераций')
            plt.grid(True, which='both', linestyle='--', alpha=0.7)
            plt.gca().invert_xaxis()
            plt.tight_layout()
            out_png = f"{prefix}_delta_iter.png"
            plt.savefig(out_png, dpi=200)
            print(f"Сохранён график: {out_png}")
            plt.show()
        else:
            print(f"Нет данных в {delta_file}")
