import matplotlib.pyplot as plt

arr1 = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001, 0.0000000001]
arr2 = [3, 6, 9, 13, 16, 19, 23, 26, 29, 33]

plt.plot(arr1, arr2)
plt.xlabel('Погрешность')
plt.ylabel('Количество итераций')
plt.grid()
plt.show()