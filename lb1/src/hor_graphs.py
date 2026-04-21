import matplotlib.pyplot as plt

if __name__ == "__main__":
    iters = [1, 1, 3, 4, 4, 5]
    eps = [1/10**i for i in range(1, 7)]
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    plt.semilogx(eps, iters, 'bo-', linewidth=2, markersize=10, markerfacecolor='red', markeredgecolor='darkred')    
    plt.xlabel('Точность Eps (логарифмическая шкала)', fontsize=14)
    plt.ylabel('Количество итераций', fontsize=14)
    plt.title('Зависимость числа итераций от точности (метод хорд)', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.7, which='both')
    
    plt.savefig("CompMath/lb1/hord.png", dpi=300)  
