import matplotlib.pyplot as plt
import numpy as np


Eps = [1/10**i for i in range(1, 7)]
Iterations = [2, 4, 7, 11, 16, 20]

plt.figure(figsize=(10, 6))
plt.semilogx(Eps, Iterations, 'bo-', linewidth=2, markersize=8)
plt.grid(True, alpha=0.3)
plt.xlabel('Точность Eps (логарифмическая шкала)')
plt.ylabel('Число итераций')
plt.title('Зависимость числа итераций от требуемой точности')
plt.axhline(y=np.log2(0.5/1e-6), color='r', linestyle='--', 
            label=f'Теоретическое значение для Eps=1e-6: {np.log2(0.5/1e-6):.1f}')
plt.legend()
plt.savefig("CompMath/lb1/bis.png", dpi=300)

