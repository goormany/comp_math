import numpy as np
import time
import random
from typing import List, Tuple

TOLERANCE = 1e-12


def get_euclidean_norm(vec: List[float]) -> float:
    """Вычисляет евклидову норму вектора"""
    sum_sq = sum(v * v for v in vec)
    return np.sqrt(sum_sq)


def get_max_norm(vec: List[float]) -> float:
    """Вычисляет максимальную норму вектора"""
    return max(abs(v) for v in vec)


def evaluate_error(mat: List[List[float]], rhs: List[float], x_curr: List[float]) -> Tuple[float, float]:
    """
    Оценивает ошибку решения (невязку)
    
    Args:
        mat: матрица коэффициентов A
        rhs: вектор правой части b
        x_curr: текущее решение
    
    Returns:
        norm_l2: L2-норма невязки
        norm_linf: L∞-норма невязки
    """
    size = len(mat)
    residual = []
    
    for row in range(size):
        ax = sum(mat[row][col] * x_curr[col] for col in range(size))
        residual.append(ax - rhs[row])
    
    norm_l2 = get_euclidean_norm(residual)
    norm_linf = get_max_norm(residual)
    
    return norm_l2, norm_linf


def init_mock_system(dim: int) -> Tuple[List[List[float]], List[float], List[float]]:
    """
    Создает тестовую систему уравнений с диагональным преобладанием
    
        dim: размерность системы
        
        mat: матрица коэффициентов A
        rhs: вектор правой части b
        x_true: точное решение
    """
    # Инициализация матрицы нулями
    mat = [[0.0] * dim for _ in range(dim)]
    rhs = [0.0] * dim
    x_true = [0.0] * dim
    
    # Генерация случайных чисел
    for i in range(dim):
        sum_off_diag = 0.0
        
        # Заполняем строку матрицы
        for j in range(dim):
            mat[i][j] = random.uniform(-15.0, 15.0)
            if i != j:
                sum_off_diag += abs(mat[i][j])
        
        # Обеспечиваем диагональное преобладание
        mat[i][i] = sum_off_diag + abs(mat[i][i]) + 2.5
        
        # Генерируем точное решение
        x_true[i] = random.uniform(-15.0, 15.0)
    
    # Вычисляем правую часть: b = A * x_true
    for i in range(dim):
        b_val = sum(mat[i][j] * x_true[j] for j in range(dim))
        rhs[i] = b_val
    
    return mat, rhs, x_true


def run_jacobi_solver(mat: List[List[float]], rhs: List[float], 
                     iteration_limit: int = 10000, required_precision: float = TOLERANCE) -> Tuple[List[float], float, float, float, bool]:
    """
    Решает СЛАУ методом Якоби (методом простых итераций)
    
    Args:
        mat: матрица коэффициентов A
        rhs: вектор правой части b
        iteration_limit: максимальное число итераций
        required_precision: требуемая точность
    
    Returns:
        solution: найденное решение
        time_spent: затраченное время
        err_l2: L2-норма невязки
        err_linf: L∞-норма невязки
        success: успешность решения
    """
    dim = len(mat)
    
    # Проверка диагональных элементов
    for i in range(dim):
        if abs(mat[i][i]) < 1e-15:
            print("[ОШИБКА] На главной диагонали найден ноль. Метод Якоби бессилен.")
            return [], 0.0, 0.0, 0.0, False
    
    # Инициализация начального приближения
    solution = [0.0] * dim
    next_solution = [0.0] * dim
    step_diff = [0.0] * dim
    
    start_time = time.time()
    
    # Итерационный процесс
    for step in range(1, iteration_limit + 1):
        # Одна итерация Якоби
        for row in range(dim):
            s = 0.0
            for col in range(dim):
                if col != row:
                    s += mat[row][col] * solution[col]
            next_solution[row] = (rhs[row] - s) / mat[row][row]
            step_diff[row] = next_solution[row] - solution[row]
        
        # Проверка сходимости
        current_error = get_euclidean_norm(step_diff)
        solution = next_solution.copy()
        
        if current_error < required_precision:
            end_time = time.time()
            time_spent = end_time - start_time
            
            err_l2, err_linf = evaluate_error(mat, rhs, solution)
            print(f"[УСПЕХ] Достигнута сходимость на {step}-й итерации.")
            return solution, time_spent, err_l2, err_linf, True
    
    # Если достигнут лимит итераций
    end_time = time.time()
    time_spent = end_time - start_time
    err_l2, err_linf = evaluate_error(mat, rhs, solution)
    
    print(f"[ПРЕДУПРЕЖДЕНИЕ] Лимит итераций ({iteration_limit}) исчерпан. Сходимость не достигнута.")
    return solution, time_spent, err_l2, err_linf, False


def main():
    size = 1000  # Уменьшил для демонстрации, можно изменить на 20000
    
    print("       РЕШЕНИЕ СЛУ МЕТОДОМ ПРОСТЫХ ИТЕРАЦИЙ      \n")
    
    print(f">> Шаг 1: Инициализация тестовой системы ({size}x{size})...")
    coeff_matrix, free_terms, exact_ans = init_mock_system(size)
    
    print(">> Шаг 2: Запуск вычислительного ядра...")
    calc_ans, duration, res_l2, res_linf, is_solved = run_jacobi_solver(coeff_matrix, free_terms)
    
    if is_solved:
        print("\n------------- СВОДКА РЕЗУЛЬТАТОВ -------------")
        print(f"Затраченное время = {duration:.6f} с")
        print(f"L_2 = {res_l2:.6e}")
        print(f"L_inf = {res_linf:.6e}")
        
        # Вычисление погрешности относительно точного решения
        absolute_error = [calc_ans[i] - exact_ans[i] for i in range(size)]
        print(f"Погрешность относительно точного решения ||x - x*||_2 = {get_euclidean_norm(absolute_error):.6e}\n")
        
        print("Детализация вектора решения (первые 5 элементов):\n")
        print(f"{'Индекс':<14} | {'Вычислено':<24} | {'Ожидалось':<30}")
        print("-" * 46)
        for i in range(min(5, size)):
            print(f"x[{i}]     | {calc_ans[i]:<15.6f} | {exact_ans[i]:<15.6f}")
        print("----------------------------------------------")
    else:
        print("\n[!] Вычисления прерваны. Проверьте матрицу на вырожденность или диагональное преобладание.")


if __name__ == "__main__":
    # Для воспроизводимости результатов
    random.seed(42)
    main()