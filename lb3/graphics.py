import matplotlib.pyplot as plt
import numpy as np

rect = []
trap = []
simp = []
epsilons = []
flag = 0
current_epsilon_index = -1  # отслеживаем текущий ε

f = open("out.txt", "r")
for line in f:
    line = line.strip()
    
    if "e = " in line:
        epsilons.append(float(line[4:]))
        current_epsilon_index += 1
        # Добавляем новые списки для каждого ε
        rect.append([])
        trap.append([])
        simp.append([])
        continue
    
    if "--" in line:
        continue
    
    if "Rectangle Method" in line:
        flag = 1
        continue
    elif "Trapezoidal Method" in line:
        flag = 2
        continue
    elif "Simpson Method" in line:
        flag = 3
        continue

    if line:  # если строка не пустая
        if flag == 1:
            rect[current_epsilon_index].append(float(line))
        elif flag == 2:
            trap[current_epsilon_index].append(float(line))
        elif flag == 3:
            simp[current_epsilon_index].append(float(line))

f.close()

rect_iterations = [len(r) for r in rect]
trap_iterations = [len(t) for t in trap]
simp_iterations = [len(s) for s in simp]

data = {
    'Rectangle Method': rect_iterations,
    'Trapezoidal Method': trap_iterations,
    'Simpson Method': simp_iterations
}

plt.figure(figsize=(10, 6))

colors = ['red', 'blue', 'green']
markers = ['o', 's', '^']

for (method, iterations), color, marker in zip(data.items(), colors, markers):
    plt.plot(epsilons, iterations, marker=marker, color=color, 
             label=method, linewidth=2, markersize=8)

plt.title('Зависимость количества итераций от ε')
plt.xlabel('ε')
plt.ylabel('Количество итераций')
plt.xscale('log')
plt.legend()
plt.grid(True, alpha=0.3)

# Добавляем подписи значений
for method, iterations in data.items():
    for eps, it in zip(epsilons, iterations):
        plt.annotate(str(it), (eps, it), textcoords="offset points", 
                    xytext=(0, 10), ha='center', fontsize=9)

plt.tight_layout()
plt.savefig("graphics.jpg")
plt.show()