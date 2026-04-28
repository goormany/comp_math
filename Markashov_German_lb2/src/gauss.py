import numpy as np
import time
from typing import List, Tuple

EPS = 1e-16

def gauss(matrix: List[List[float]]) -> Tuple[List[float], bool]:
    """
    Метод Гаусса с выбором главного элемента по столбцу
    
    Args:
        matrix: расширенная матрица [A|b] размера n x (n+1)
    
    Returns:
        x: вектор решений
        success: True если решение найдено, False если матрица вырождена
    """
    # Создаем копию матрицы, чтобы не изменять оригинал
    A = [row[:] for row in matrix]
    n = len(A)
    x = [0.0] * n
    
    # Прямой ход
    for i in range(n):
        # Поиск главного элемента
        max_row = i
        max_val = abs(A[i][i])
        
        for j in range(i + 1, n):
            if abs(A[j][i]) > max_val:
                max_val = abs(A[j][i])
                max_row = j
        
        # Проверка на вырожденность
        if max_val < EPS:
            return x, False
        
        # Перестановка строк
        if max_row != i:
            A[i], A[max_row] = A[max_row], A[i]
        
        # Исключение переменных
        pivot = A[i][i]
        for j in range(i + 1, n):
            factor = A[j][i] / pivot
            for k in range(i, n + 1):
                A[j][k] -= factor * A[i][k]
    
    # Обратный ход
    for i in range(n - 1, -1, -1):
        sum_ax = 0.0
        for j in range(i + 1, n):
            sum_ax += A[i][j] * x[j]
        x[i] = (A[i][n] - sum_ax) / A[i][i]
    
    return x, True


def check_error(A: List[List[float]], b: List[float], x: List[float]) -> Tuple[float, float]:
    """
    Вычисляет невязку решения
    
    Args:
        A: матрица коэффициентов
        b: вектор правой части
        x: найденное решение
    
    Returns:
        l2_err: L2-норма невязки
        linf_err: L∞-норма невязки
    """
    n = len(A)
    l2_err = 0.0
    linf_err = 0.0
    
    for i in range(n):
        # Вычисляем сумму A[i][j] * x[j]
        sum_ax = sum(A[i][j] * x[j] for j in range(n))
        
        resid = sum_ax - b[i]
        
        l2_err += resid * resid
        linf_err = max(linf_err, abs(resid))
    
    l2_err = np.sqrt(l2_err)
    return l2_err, linf_err


def run_test(matrix: List[List[float]], test_name: str, print_details: bool = False):
    """
    Запускает тест для заданной матрицы
    
    Args:
        matrix: расширенная матрица [A|b]
        test_name: название теста
        print_details: печатать ли подробные результаты
    """
    n = len(matrix)
    
    # Разделяем матрицу A и вектор b
    A = [row[:-1] for row in matrix]
    b = [row[-1] for row in matrix]
    
    # Решение СЛАУ
    start_time = time.time()
    x, success = gauss(matrix)
    end_time = time.time()
    
    comp_time = end_time - start_time
    
    print(f"=== Результаты: {test_name} ===")
    
    if not success:
        print("ОШИБКА: Матрица вырождена.\n")
        return
    
    # Вычисляем невязки
    norm_l2, norm_linf = check_error(A, b, x)
    
    if print_details:
        for i in range(n):
            print(f"x[{i + 1}] = {x[i]:.5f}")
        
        print("\nПодробная обратная проверка (A * x = b):")
        for i in range(n):
            calc_b = sum(A[i][j] * x[j] for j in range(n))
            diff = abs(calc_b - b[i])
            print(f"Уравнение {i + 1}:")
            terms = []
            for j in range(n):
                terms.append(f"({A[i][j]} * {x[j]:.5f})")
            print("  " + " + ".join(terms))
            print(f"  = {calc_b:.10f}  (Ожидалось: {b[i]}, Разница: {diff:.2e})\n")
    else:
        print("Вывод корней сокращен (матрица слишком велика):")
        print(f"x[1] = {x[0]:.5f}")
        print(f"x[2] = {x[1]:.5f}")
        print("...")
        print(f"x[{n}] = {x[-1]:.5f}")
        print("(Подробная обратная проверка пропущена для экономии времени)\n")
    
    print(f"Невязка L2 (Евклидова) : {norm_l2:.2e}")
    print(f"Невязка L-inf (Макс)   : {norm_linf:.2e}")
    print(f"Время работы          : {comp_time:.6f} сек.\n")


def run_file_test(filename: str):
    """
    Читает матрицу из файла и запускает тест
    
    Args:
        filename: имя файла с матрицей
    """
    try:
        with open(filename, 'r') as fin:
            # Читаем размер матрицы
            n = int(fin.readline().strip())
            
            # Читаем матрицу
            matrix = []
            for _ in range(n):
                row = list(map(float, fin.readline().split()))
                if len(row) != n + 1:
                    print(f"ОШИБКА: Неверный формат строки. Ожидалось {n + 1} чисел")
                    return
                matrix.append(row)
        
        # Решаем систему
        print_details = False  # (n <= 20) - можно изменить при необходимости
        run_test(matrix, f"Матрица {n}x{n} из файла", print_details)
        
    except FileNotFoundError:
        print(f"=== Результаты: Тест из файла {filename} ===")
        print(f"ОШИБКА: Не удалось открыть файл {filename}\n")
    except ValueError as e:
        print(f"ОШИБКА: Неверный формат данных в файле - {e}\n")


def main():
    # 1. Тест матрицы Гильберта
    n_hilbert = 10
    hilbert = []
    for i in range(n_hilbert):
        row = []
        # Заполняем матрицу Гильберта
        for j in range(n_hilbert):
            row.append(1.0 / (i + j + 1))
        # Сумма строки для правой части
        row_sum = sum(row)
        row.append(row_sum)
        hilbert.append(row)
    
    run_test(hilbert, "Матрица Гильберта 10x10", print_details=True)
    
    # 2. Тест хорошо обусловленной СЛУ
    simple_sys = [
        [4.0, -1.0, 0.0, 3.0],
        [-1.0, 4.0, -1.0, 2.0],
        [0.0, -1.0, 4.0, 3.0]
    ]
    run_test(simple_sys, "Хорошо обусловленная СЛУ 3x3", print_details=True)
    
    # 3. Тест матрицы из файла
    run_file_test("input.txt")


if __name__ == "__main__":
    main()