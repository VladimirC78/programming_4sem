import numpy as np
import matplotlib.pyplot as plt
import random as rd
from mpl_toolkits.mplot3d import Axes3D
import time

animation = False

def create_distributed_points(center, matrix, d):
    """
    Генерирует случайную точки внутри заданного эллипсоида.
    
    Параметры:
        center (np.array): Центр эллипсоида [x, y, z]
        matrix (np.array): Матрица 3x3, задающая форму эллипсоида
        d (int): размерность пространства
        
    Возвращает:
        np.array: точку в форме (d, 1)
    """

    # Полуоси и их длина

    evalues, half_axes = np.linalg.eigh(matrix)
    radii = evalues**(-0.5)

    ratio = np.random.uniform(0, 1)

    angles = np.array([])
    point = np.array([])
    
    # Создаем массив углов

    for i in range(d - 1):
        if i == 0:
            angles = np.append(angles, rd.uniform(0, 2*np.pi))
        else:
            angles = np.append(angles, rd.uniform(0, np.pi))
    
    # Переводим координаты точки из гиперсферических в декартовы

    for i in range(d - 1):
        x = radii[i] * ratio
        for j in range(d - i - 1):
            if i == 0:
                x *= np.sin(angles[j])
            else:
                x *= np.cos(angles[i - 1])
                if i + j < d - 1:
                    x *= np.sin(angles[i + j])
        point = np.append(point, x)
    
    point = np.append(point, radii[-1] * np.cos(angles[-1]))

    # И поворачиваем в соответствии с направлениями полуосей

    point = half_axes @ point
    return center + point


def locate_farthest(X, n):
    """
    Возвращает наиболее удаленную точку в пространстве R^{d + 1}.
    
    Параметры:
        X (np.array): Матрица координат в R^{d + 1}
        n (int): Число точек
        
    Возвращает:
        np.array: Наиболее удаленную точку в форме (d + 1, 1)
    """
        
    ind = 0
    max_norm = 0

    for i in range(n):
        norm = np.linalg.norm(X.T[i])
        if norm > max_norm:
            max_norm = norm
            ind = i

    return X.T[ind]


def build_contraction_operator(xi, alpha, d):
    """
    Возвращает оператор сжатия в R^{d + 1} вдоль заданного вектора
    
    Параметры:
        xi (np.array): Вектор, задающий направление сжатия
        alpha (float): Коэффициент сжатия
        d (int): Размерность пространства
        
    Возвращает:
        np.array: Матрицу оператора сжатия в форме (d + 1, d + 1)
    """

    uni_vec = xi / np.linalg.norm(xi)

    xixi = np.dot(uni_vec.reshape(d + 1, 1), uni_vec.reshape(1, d + 1))

    R = np.eye(d + 1) - (1 - alpha)* xixi
    return R


def form_coord_matrix(n, d):
    """
    Возвращает матрицу координат точек в гиперплоскости x_{d + 1} = 1.
    
    Параметры:
        n (int): Число точек
        d (int): Размерность пространства
        
    Возвращает:
        np.array: Координатную матрицу формы (n, d + 1)
    """
    X = np.array([])
    det = 0

    # Создаем симметричную положительно определенную матрицу эллипса, в котором раскиданы точки
    
    A = np.eye(d) - np.eye(d)
    while det == 0:
        A = (np.random.rand(d, d) - 0.5) * 0.15
        det = np.linalg.det(A)
    test_matrix = A.T @ A 

    center = np.array([0] * d)

    # Раскидываем точки внутри полученного выше эллипса

    for i in range(n):
        point = np.append(create_distributed_points(center, test_matrix, d), 1)
        if i == 0:  X = point
        else:   X = np.vstack((X, point))


    return X.T, test_matrix

def form_main_minor(A_m):
    """
    Выделяет главный минор матрицы.
    
    Параметры:
        A_m (np.array): Матрица
        
    Возвращает:
        np.array: Главный минор матрицы
    """
    B = (A_m[:-1].T)[:-1]
    return B

def build_ellipse(center, ell_matrix, num_points = 1000):
    """
    Создает данные для отрисовки двумерного эллипса по его матрице и центру

    Параметры:
        center (np.array): Вектор координат центра эллипса
        ell_matrix (np.array): Матрица эллипса формы (2, 2)

    Возвращает:
        np.array: Данные для построения эллипса формы (2, num_points)
    """
    evalues, half_axes = np.linalg.eigh(ell_matrix)
    radii = evalues**(-0.5)
    theta = np.linspace(0, 2*np.pi, num_points)
    ellipse_set = np.array([radii[0] * np.cos(theta), radii[1] * np.sin(theta)])
    trans_ellipse = half_axes @ ellipse_set
    output = np.array([trans_ellipse[0] + center[0], trans_ellipse[1] + center[1]])
    return output


def build_ellipsoid(center, matrix, num_points = 1000):
    """
    Создает данные для отрисовки трехмерного эллипсоида по его матрице и центру

    Параметры:
        center (np.array): Вектор координат центра эллипса
        ell_matrix (np.array): Матрица эллипсоида формы (3, 3) 
    Возвращает:
        np.array: Данные для построения эллипса формы (3, num_points)
    """

    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    radii = 1 / np.sqrt(eigenvalues)

    u = np.linspace(0, 2 * np.pi, num_points)
    v = np.linspace(0, np.pi, num_points)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))

    points = np.stack([x * radii[0], y * radii[1], z * radii[2]])
    points_rotated = np.dot(eigenvectors, points.reshape(3, -1)).reshape(3, num_points, num_points)

    x_final = points_rotated[0] + center[0]
    y_final = points_rotated[1] + center[1]
    z_final = points_rotated[2] + center[2]

    output = np.array([x_final, y_final, z_final])

    return output

def perform_step(step, X, cont_seq, n, d):
    """
    Осуществляет одну итерацию алгоритма
    """

    # Коэффициент сжатия
    beta = 1 / (step + 2)
    alpha = 1 - beta

    # Находим самую удаленную точку 
    farthest = locate_farthest(X, n)

    # Создаем оператор сжатия по направлению этой точки
    R = build_contraction_operator(farthest, alpha, d)

    # Расчитываем новые координаты точек
    X = R @ X

    # Обновляем оператор перехода от изначального пространства
    cont_seq = R @ cont_seq

    return X, cont_seq


def find_ellipse(cont_seq, X, n):
    """
    Находит матрицу и центр эллипса как сечения вспомогательного эллипсоида

    Параметры:
        cont_seq (np.array): матрица преобразования от исходного пространства
        X (np.array): матрица координат точек в сжатом пространстве
        n (int): количество точек
    
    Возвращает:
        np.array, np.array: вектор координат центра и матрицу минимального эллипса
    """

    # Матрица вспомогательного минимального эллипсоида в R^{d + 1}

    A_m = (cont_seq.T @ cont_seq) / np.linalg.norm(locate_farthest(X, n))**2

    # Берем его сечение гиперплоскостью x_{d + 1} = 1

    B = form_main_minor(A_m)

    b = A_m.T[-1][:-1]

    center = -np.dot(np.linalg.inv(B), b)

    ell_matrix = B / (1 - np.dot(center, b) - A_m[-1][-1])

    return center, ell_matrix


def draw2d (center, ell_matrix, test_matrix, A0, step = 0):
    plt.figure(figsize=(6, 6))
    plt.errorbar(A0[0], A0[1], fmt=".k")
    ellipse = build_ellipse(center, ell_matrix)
    test_ell = build_ellipse([0, 0], test_matrix)
    if not animation:
        plt.plot(ellipse[0], ellipse[1], "r", label="Approximate solution, step " + str(m))
    else:
        plt.plot(ellipse[0], ellipse[1], "r", label="Approximate solution, step " + str(step))
    plt.plot(test_ell[0], test_ell[1], "b", label="Initial ellipse")
    plt.legend()

    if animation:
        plt.savefig('general_problem/animation/2Dstep' + (str(step) if step > 9 else '0' + str(step)) + '.png')

def draw3d (center, ell_matrix, A0, step = 0):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ellipsoid = build_ellipsoid(center, ell_matrix, 1000)

    ax.errorbar(A0[0], A0[1], A0[2], xerr = 0, yerr = 0, zerr = 0, fmt='.k')
    ax.plot_surface(ellipsoid[0], ellipsoid[1], ellipsoid[2], color="r", alpha=0.3)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    if animation:
        plt.savefig('general_problem/animation/3Dstep' + (str(step) if step > 9 else '0' + str(step)) + '.png')


print("Enter number of dimensions")

d = int(input())

# print("Enter number of points")

# n = int(input())

# print("Enter number of steps")

# m = int(input())

m = 100

f = open('time_measurements_n.txt', 'w')

for n in range(10, 100000, 100):

    start_time = time.time()

    X, test_matrix = form_coord_matrix(n, d)

    A0 = X

    cont_seq = np.eye(d + 1)

    # Решаем задачу в пространстве R^{d + 1} для эллипсоида с фиксированным центром

    for step in range(m):
        X, cont_seq = perform_step(step, X, cont_seq, n, d)
        if animation:
            center, ell_matrix = find_ellipse(cont_seq, X, n)
            if d == 2:
                draw2d(center, ell_matrix, test_matrix, A0, step)
            elif d == 3:
                draw3d(center, ell_matrix, A0, step)
    
    end_time = time.time()
    t = end_time - start_time

    f.write(str(n) + ' ' + str(t) + '\n')
    print(n)

f.close()



if not animation:
    center, ell_matrix = find_ellipse(cont_seq, X, n)

    # Нарисуем это все для случаев размерности 2 или 3

    if d == 2:
        draw2d(center, ell_matrix, test_matrix, A0)
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