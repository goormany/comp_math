import matplotlib.pyplot as plt

if __name__ == "__main__":
    iters = [4, 7, 10, 14, 17, 20]
    eps = [1/10**i for i in range(1, 7)]
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    plt.semilogx(eps, iters, 'bo-', linewidth=2, markersize=10, markerfacecolor='red', markeredgecolor='darkred')    
    plt.xlabel('Точность Eps (логарифмическая шкала)', fontsize=14)
    plt.ylabel('Количество итераций', fontsize=14)
    plt.title('Зависимость числа итераций от точности (метод бисекции)', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.7, which='both')
    
    plt.savefig("bisec.png", dpi=300)  
    
    
    fig2, ax2 = plt.subplots(figsize=(12, 8))
    deltas = [1/10**i for i in range(1, 7)]
    roots = [-0.75, -0.75, -0.754883, -0.754883, -0.754883, -0.754868]
    true_root = -0.754877
    deviations = [abs(r - true_root) for r in roots]
    
    ax2.semilogx(deltas, deviations, 'ro-', linewidth=2, markersize=12,
                markerfacecolor='red', markeredgecolor='darkred')
    
    ax2.set_xlabel('Delta', fontsize=14)
    ax2.set_ylabel('Отклонение от точного корня', fontsize=14)
    ax2.set_title('Зависимость точности корня от ошибок в исходных данных', fontsize=16)
    ax2.grid(True, linestyle='--', alpha=0.7, which='both')
    plt.savefig("bisec_deviation.png", dpi=300)
    