import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("interpolation_data.csv")
nodes = pd.read_csv("nodes.csv")

x = data["x"]
fx = data["f(x)"]
lagrange = data["lagrange"]
spline = data["spline"]

xn = nodes["x"]
yn = nodes["y"]

plt.figure(figsize=(12, 7))

plt.plot(x, fx, linewidth=2, label="Исходная функция f(x)")
plt.plot(x, lagrange, linewidth=2, label="Интерполяционный многочлен Лагранжа")
plt.plot(x, spline, linewidth=2, label="Закреплённый кубический сплайн")

plt.scatter(
    xn, yn,
    s=30,
    marker="o",
    label="Узлы интерполяции",
    zorder=5
)

plt.title("Интерполяция функции arctg(x) - 1/x на отрезке [0.1; 5]")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig("interpolation_plot.png", dpi=200)
plt.show()