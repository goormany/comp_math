import matplotlib.pyplot as plt

if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(12, 8))
    iters = [0, 3, 7, 10, 13, 17, 20, 23, 27, 30]
    iters_hord = [3, 7, 10, 13, 16, 20, 23, 26, 30, 33]
    iters_newton = [2, 3, 4, 5, 5, 5, 5, 6, 6, 6]
    iters_iter = [2, 2, 2, 2, 2, 9, 12, 14, 17, 19]
    eps = [1/10**i for i in range(1, 11)]
    
    
    
    plt.semilogx(eps, iters_iter, 'bo-', linewidth=2, markersize=10, markerfacecolor='red', markeredgecolor='darkred')    
    plt.xlabel('Точность Eps (логарифмическая шкала)', fontsize=14)
    plt.ylabel('Количество итераций', fontsize=14)
    plt.title('Зависимость числа итераций от точности (метод Ньютона)', fontsize=16)
    plt.title('Зависимость числа итераций от точности (метод Простых итераций)', fontsize=16)
    # plt.title('Зависимость числа итераций от точности (метод бисекции)', fontsize=16)
    # plt.title('Зависимость числа итераций от точности (метод хорд)', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.7, which='both')
    # plt.show()
    # plt.savefig("hord.png", dpi=300)
    # plt.savefig("bisec.png", dpi=300)
    # plt.savefig("Newton.png", dpi=300)
    plt.savefig("iter.png", dpi=300)