import math

def f(x):
    return math.cos(x) * math.exp(-x**2)

def rectangle_method(a, b, n):
    """Метод средних прямоугольников"""
    h = (b - a) / n
    total_sum = 0.0
    for i in range(n):
        # Берем значение в середине каждого отрезка
        total_sum += f(a + h * (i + 0.5))
    return h * total_sum

def trapezoid_method(a, b, n):
    """Метод трапеций"""
    h = (b - a) / n
    total_sum = (f(a) + f(b)) / 2.0
    for i in range(1, n):
        total_sum += f(a + i * h)
    return h * total_sum

def simpson_method(a, b, n):
    """Метод Симпсона (парабол)"""
    if n % 2 != 0: n += 1  # n должно быть четным
    h = (b - a) / n
    total_sum = f(a) + f(b)
    
    for i in range(1, n):
        x = a + i * h
        # Коэффициенты 4, 2, 4, 2...
        if i % 2 == 0:
            total_sum += 2.0 * f(x)
        else:
            total_sum += 4.0 * f(x)
            
    return h * total_sum / 3.0

def runge_refinement(name, method, a, b, eps, order):
    """
    Автоматический подбор шага по правилу Рунге.
    order: порядок точности метода (2 для прямоуг./трапеций, 4 для Симпсона).
    """
    n = 2
    # Знаменатель в формуле Рунге: 2^p - 1
    k = pow(2, order) - 1
    
    while True:
        i1 = method(a, b, n)
        i2 = method(a, b, 2 * n)
        
        # Оценка погрешности
        error = abs(i2 - i1) / k
        # Уточненный результат по Рунге
        result = i2 + (i2 - i1) / k
        
        if error <= eps:
            break
        n *= 2
        
    print(f"{name}")
    print(f"I = {result:.10f}")
    print(f"n = {2 * n}")  # итоговое количество разбиений
    print(f"h = {(b - a) / (2 * n):.10f}")
    print(f"Погрешность по Рунге = {error:.10e}")
    print("-" * 30)

def main():
    print("Вариант 9")
    print("Интеграл: cos(x)*exp(-x^2) на [0, 1]")
    
    try:
        eps = float(input("Введите точность eps (например, 0.0001): "))
    except ValueError:
        print("Ошибка ввода. Используйте число.")
        return

    a, b = 0, 1
    
    runge_refinement("Метод прямоугольников:", rectangle_method, a, b, eps, 2)
    runge_refinement("Метод трапеций:", trapezoid_method, a, b, eps, 2)
    runge_refinement("Метод Симпсона:", simpson_method, a, b, eps, 4)

if __name__ == "__main__":
    main()