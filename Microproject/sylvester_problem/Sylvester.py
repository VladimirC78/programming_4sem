import numpy as np
import cvxopt as cv
import random as rd
import matplotlib.pyplot as plt

n = 1000

A0 = []

desired = 20

for i in range(n):
    x = 1000
    y = 1000
    while (x ** 2 + y ** 2 > desired ** 2):
        x = rd.uniform(-20, 20)
        y = rd.uniform(-20, 20)
    A0.append([x, y])

A0 = np.array(A0)

b0 = []

for i in range(n):
    sum = 0
    for j in range(len(A0[i])):
        sum += 0.5 * A0[i][j]**2
    b0.append(sum)


p = cv.matrix(-np.array(b0))

Q = cv.matrix(np.dot(A0, A0.T))

G = cv.matrix(-np.eye(n))

h = cv.matrix(np.zeros(n))

A = cv.matrix(np.array([1.0] * n), (1, n))

b = cv.matrix(1.0)


solution = cv.solvers.qp(Q, p, G, h, A, b)

u = solution['x']

u = u.T

center = np.dot(u, A0)


f = 0.5 * np.linalg.norm(np.dot(u, A0))**2 + u * p

radius = (-2 * f)**0.5

theta = np.linspace(0, 2 * np.pi, 1000)

x = radius * np.cos(theta)
y = radius * np.sin(theta)

print(radius)

print(center)

x += center[0][0]
y += center[0][1]


plt.figure(figsize=(6, 6))
plt.errorbar(A0.T[0], A0.T[1], fmt=".k")
plt.plot(x[0], y[0], "r", linewidth=2)
plt.grid()
plt.show()