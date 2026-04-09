import matplotlib.pyplot as plt
import numpy as np


x_data = np.linspace(-5, 5, 100)
x_data = x_data[x_data != 0]

y_data = np.arctan(1 / x_data) - x_data ** 2

plt.figure(figsize=(10, 6))
plt.plot(x_data, y_data, linewidth=2)
plt.grid(True, alpha=0.3)
plt.xticks(np.arange(-5, 6, 1)) 
plt.xlabel("X")
plt.ylabel("Y")
plt.title(r"График функции: $\arctan(1/x) - x^2$")
plt.axhline(0, color="black", linewidth=0.5)
plt.axvline(0, color="black", linewidth=0.5)
plt.show()

# x_data = np.linspace(-5, 5, 100)
# y_data = - (1 / (x_data**2 + 1)) - 2*x_data

# plt.figure(figsize=(10, 6))
# plt.plot(x_data, y_data, linewidth=2)
# plt.grid(True, alpha=0.3)
# plt.xticks(np.arange(-5, 6, 1)) 
# plt.xlabel("X")
# plt.ylabel("Y")
# plt.title(r"График функции: $-\frac{1}{x^2+1} - 2x$")
# plt.axhline(0, color="black", linewidth=0.5)
# plt.axvline(0, color="black", linewidth=0.5)
# plt.show()
# Для f(x) = arctan(1/x) - x**2:
# x != 0 (разрыв)
# Функция непрерывна на (-беск, 0) и (0, +беск)

#  f'(x) = -1/(x²(1+x²)) - 2x

eps = np.array([0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001])
N_bysection = np.array([0, 3, 6, 9, 13, 16])
N_chord = np.array([2, 3, 4, 5, 6, 7])
N_newton = np.array([1, 2, 2, 3, 3, 4])
N_iter = np.array([1, 2, 3, 3, 4, 4])

plt.figure(figsize=(8, 5))
plt.plot(eps, N_bysection, 'bo-', label='bysection', linewidth=2)
plt.plot(eps, N_chord, 'rs-', label='horda', linewidth=2)
plt.plot(eps, N_newton, 'g^-', label='newton', linewidth=2)
plt.plot(eps, N_newton, 'y*-', label='iter', linewidth=2)


plt.grid(True, alpha=0.3)
plt.xlabel("Точность Eps")
plt.ylabel("Число итераций N")
plt.title("Сравнение сходимости методов")
plt.show()
