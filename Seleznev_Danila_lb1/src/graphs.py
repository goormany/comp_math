import matplotlib.pyplot as plt
import pandas as pd

# Чтение данных из CSV файлов
convergence_df = pd.read_csv('convergence_data.csv')
sensitivity_df = pd.read_csv('sensitivity_data.csv')

# Извлечение данных для каждого метода
bisect_conv = convergence_df[convergence_df['Method'] == 'Bisect']
horda_conv = convergence_df[convergence_df['Method'] == 'Horda']
newton_conv = convergence_df[convergence_df['Method'] == 'Newton']
iter_conv = convergence_df[convergence_df['Method'] == 'Iter']

bisect_sens = sensitivity_df[sensitivity_df['Method'] == 'Bisect']
horda_sens = sensitivity_df[sensitivity_df['Method'] == 'Horda']
newton_sens = sensitivity_df[sensitivity_df['Method'] == 'Newton']
iter_sens = sensitivity_df[sensitivity_df['Method'] == 'Iter']

# =========================
# График 1: Бисекция - скорость сходимости
# =========================
plt.figure(figsize=(8,5))
plt.plot(bisect_conv['Eps'], bisect_conv['Iterations'], marker='o', linestyle='-', color='blue', label='Бисекция')
plt.xscale('log')
plt.xlabel('Eps (точность)')
plt.ylabel('Число итераций')
plt.title('Метод бисекции: зависимость числа итераций от Eps')
plt.grid(True)
plt.legend()
plt.show()

# =========================
# График 2: Бисекция - устойчивость к ошибкам
# =========================
true_root = 0.154406

plt.figure(figsize=(8,5))
plt.plot(bisect_sens['Delta'], abs(bisect_sens['Root'] - true_root), marker='o', linestyle='-', color='blue', label='Бисекция')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta (ошибка округления)')
plt.ylabel('|Корень - истинное значение|')
plt.title('Метод бисекции: устойчивость к ошибкам')
plt.grid(True)
plt.legend()
plt.show()

# =========================
# График 3: Хорда - скорость сходимости
# =========================
plt.figure(figsize=(8,5))
plt.plot(horda_conv['Eps'], horda_conv['Iterations'], marker='s', linestyle='-', color='green', label='Хорда')
plt.xscale('log')
plt.xlabel('Eps (точность)')
plt.ylabel('Число итераций')
plt.title('Метод хорд: зависимость числа итераций от Eps')
plt.grid(True)
plt.legend()
plt.show()

# =========================
# График 4: Хорда - устойчивость к ошибкам
# =========================
plt.figure(figsize=(8,5))
plt.plot(horda_sens['Delta'], abs(horda_sens['Root'] - true_root), marker='s', linestyle='-', color='green', label='Хорда')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta (ошибка округления)')
plt.ylabel('|Корень - истинное значение|')
plt.title('Метод хорд: устойчивость к ошибкам')
plt.grid(True)
plt.legend()
plt.show()

# =========================
# График 5: Ньютон - скорость сходимости
# =========================
plt.figure(figsize=(8,5))
plt.plot(newton_conv['Eps'], newton_conv['Iterations'], marker='^', linestyle='-', color='red', label='Ньютон')
plt.xscale('log')
plt.xlabel('Eps (точность)')
plt.ylabel('Число итераций')
plt.title('Метод Ньютона: зависимость числа итераций от Eps')
plt.grid(True)
plt.legend()
plt.show()

# =========================
# График 6: Ньютон - устойчивость к ошибкам
# =========================
plt.figure(figsize=(8,5))
plt.plot(newton_sens['Delta'], abs(newton_sens['Root'] - true_root), marker='^', linestyle='-', color='red', label='Ньютон')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta (ошибка округления)')
plt.ylabel('|Корень - истинное значение|')
plt.title('Метод Ньютона: устойчивость к ошибкам')
plt.grid(True)
plt.legend()
plt.show()

# =========================
# График 7: Простые итерации - скорость сходимости
# =========================
plt.figure(figsize=(8,5))
plt.plot(iter_conv['Eps'], iter_conv['Iterations'], marker='d', linestyle='-', color='purple', label='Простые итерации')
plt.xscale('log')
plt.xlabel('Eps (точность)')
plt.ylabel('Число итераций')
plt.title('Метод простых итераций: зависимость числа итераций от Eps')
plt.grid(True)
plt.legend()
plt.show()

# =========================
# График 8: Простые итерации - устойчивость к ошибкам
# =========================
plt.figure(figsize=(8,5))
plt.plot(iter_sens['Delta'], abs(iter_sens['Root'] - true_root), marker='d', linestyle='-', color='purple', label='Простые итерации')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta (ошибка округления)')
plt.ylabel('|Корень - истинное значение|')
plt.title('Метод простых итераций: устойчивость к ошибкам')
plt.grid(True)
plt.legend()
plt.show()

# =========================
# График 9: Сравнение скорости сходимости всех методов
# =========================
plt.figure(figsize=(10,6))
plt.plot(bisect_conv['Eps'], bisect_conv['Iterations'], marker='o', linestyle='-', color='blue', label='Бисекция')
plt.plot(horda_conv['Eps'], horda_conv['Iterations'], marker='s', linestyle='-', color='green', label='Хорда')
plt.plot(newton_conv['Eps'], newton_conv['Iterations'], marker='^', linestyle='-', color='red', label='Ньютон')
plt.plot(iter_conv['Eps'], iter_conv['Iterations'], marker='d', linestyle='-', color='purple', label='Простые итерации')
plt.xscale('log')
plt.xlabel('Eps (точность)')
plt.ylabel('Число итераций')
plt.title('Сравнение скорости сходимости методов')
plt.grid(True)
plt.legend()
plt.show()

# =========================
# График 10: Устойчивость методов - отклонение корня от истинного значения
# =========================
true_root = 0.154406

plt.figure(figsize=(10,6))
plt.plot(bisect_sens['Delta'], abs(bisect_sens['Root'] - true_root), marker='o', linestyle='-', color='blue', label='Бисекция')
plt.plot(horda_sens['Delta'], abs(horda_sens['Root'] - true_root), marker='s', linestyle='-', color='green', label='Хорда')
plt.plot(newton_sens['Delta'], abs(newton_sens['Root'] - true_root), marker='^', linestyle='-', color='red', label='Ньютон')
plt.plot(iter_sens['Delta'], abs(iter_sens['Root'] - true_root), marker='d', linestyle='-', color='purple', label='Простые итерации')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta (ошибка округления)')
plt.ylabel('|Корень - истинное значение|')
plt.title('Устойчивость методов: отклонение корня от истинного значения')
plt.grid(True)
plt.legend()
plt.show()
