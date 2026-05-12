import matplotlib.pyplot as plt

if __name__ == "__main__":
    iters = [2, 5, 9, 13, 16, 20]
    eps = [1/10**i for i in range(1, 7)]

    plt.semilogx(eps, iters, 'bo-', linewidth=2, markersize=10, markerfacecolor='red', markeredgecolor='darkred')    
    plt.xlabel('Точность Eps (логарифмическая шкала)', fontsize=14)
    plt.ylabel('Количество итераций', fontsize=14)
    plt.title('Зависимость числа итераций от точности (метод iter)', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.7, which='both')
    
    plt.savefig("iter.png", dpi=300)  
    
    fig2, ax2 = plt.subplots(figsize=(12, 8))
    deltas = [1/10**i for i in range(1, 5)]
    roots = [-0.755, -0.7549, -0.75487, -0.754887]
    true_root = -0.754877
    deviations = [abs(r - true_root) for r in roots]
    
    ax2.semilogx(deltas, deviations, 'ro-', linewidth=2, markersize=12,
                markerfacecolor='red', markeredgecolor='darkred')
    
    ax2.set_xlabel('Delta', fontsize=14)
    ax2.set_ylabel('Отклонение от точного корня', fontsize=14)
    ax2.set_title('Зависимость точности корня от ошибок в исходных данных', fontsize=16)
    ax2.grid(True, linestyle='--', alpha=0.7, which='both')
    plt.savefig("iter_deviation.png", dpi=300)