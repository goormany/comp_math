import random

def generate_matrix(filename, n):
    print(f"Начинаю генерацию матрицы {n}x{n}...")
    
    with open(filename, 'w') as f:
        f.write(f"{n}\n")
        
        for i in range(n):
            row = []
            sum_b = 0.0 
            
            for j in range(n):
                if i == j:
                    val = random.uniform(15000.0, 20000.0)
                else:
                    val = random.uniform(-10.0, 10.0)
                
                row.append(val)
                sum_b += val
            
            row.append(sum_b)
            
            f.write(" ".join(f"{x:.4f}" for x in row) + "\n")
            
    print(f"Успешно! Матрица записана в файл: {filename}")

generate_matrix("input.txt", 1000)