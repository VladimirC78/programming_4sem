import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from functions import *


print("Enter number of dimensions")

d = int(input())

print("Enter number of points")

n = int(input())

print("Enter number of steps")

m = int(input())


X, test_matrix = form_coord_matrix(n, d)

A0 = X

cont_seq = np.eye(d + 1)

# Решаем задачу в пространстве R^{d + 1} для эллипсоида с фиксированным центром

for step in range(m):
    X, cont_seq = perform_step(step, X, cont_seq, n, d)


center, ell_matrix = find_ellipse(cont_seq, X, n)

# Нарисуем это все для случаев размерности 2 или 3

if d == 2:
    draw2d(center, ell_matrix, test_matrix, A0, m)
elif d == 3:
    draw3d(center, ell_matrix, A0)

if d == 2 or d == 3:
    plt.grid()
    plt.legend()
    plt.show()


# Или хотя бы выведем координату центра и матрицу эллипсоида, чтобы был какой-то результат для d > 3 :)

print("Координаты центра")
print(center)

print("Матрица эллипсоида")
print(ell_matrix)

print("Матрица исходного эллипсоида")
print(test_matrix)

print("Разность матриц")
print(ell_matrix - test_matrix)