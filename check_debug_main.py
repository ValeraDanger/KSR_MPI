import numpy as np
import os

# Активировать отладку (True - включить, False - выключить)
debug = True

# Создаем файл для отладочной информации
if debug:
    debug_file = open("py_debug.txt", "w")
    np.set_printoptions(precision=10, suppress=True)
    debug_file.write("=== Метод серединных градиентов (Python): отладочная информация ===\n")

A = np.array(
[ 
    [-4, 1, 1, 0],
    [1, -4, 0, 1],
    [1, 0, -4, 1],
    [0, 1, 1, -4]
]
)
A = 9*A
x0 = np.zeros(4).reshape(-1, 1)

b = [-7.535398163397445, -7.646055384571343, -0.897879165856758, -1.0085363870306536]
#b =   [-4415.278150, -17325.120336,  -880.528004, -7577.412018,  -613.908427, -5110.962206,  -301.706564,  -354.993071,  -636.136165,    22.222222, -3269.246230,  -111.825866,   -85.942352,  -149.007995,  -273.802169, -2540.699799]
b = np.array(b).reshape(-1, 1)

# Итерация 0
if debug:
    debug_file.write("Итерация 0:\n")
    debug_file.write(f"x0 = {x0.flatten().tolist()}\n")

h0 = -b.copy()

if debug:
    debug_file.write(f"h0 = {h0.flatten().tolist()}\n")
    debug_file.write(f"A @ h0 = {(A @ h0).flatten().tolist()}\n")

print(np.shape(A @ x0))

# Вычисление alpha0
alpha = - (np.dot((A @ x0 - b).T, h0)) / (np.dot((A @ h0).T, h0))
print(alpha)

if debug:
    debug_file.write(f"alpha0 = {alpha[0][0]}\n")

# Итерация 1
x1 = x0 + alpha*h0

if debug:
    debug_file.write(f"x1 = {x1.flatten().tolist()}\n")

r1 = A @ x1 - b

if debug:
    debug_file.write(f"r1 = {r1.flatten().tolist()}\n")

# Вычисление beta0
beta = (np.dot((A @ h0).T, r1)) / (np.dot((A @ h0).T, h0))

if debug:
    debug_file.write(f"beta0 = {beta[0][0]}\n")

# Обновление направления
h1 = -r1 + beta * h0

if debug:
    debug_file.write(f"h1 = {h1.flatten().tolist()}\n")
    debug_file.write("-------------------------------------\n")
    debug_file.write(f"A @ h1 = {(A @ h1).flatten().tolist()}\n")

# Вычисление alpha1
alpha = - (np.dot((A @ x1 - b).T, h1)) / (np.dot((A @ h1).T, h1))

if debug:
    debug_file.write(f"alpha1 = {alpha[0][0]}\n")

# Итерация 2
x2 = x1 + alpha * h1
print(x2)

if debug:
    debug_file.write(f"x2 = {x2.flatten().tolist()}\n")
    r2 = A @ x2 - b
    debug_file.write(f"r2 = {r2.flatten().tolist()}\n")
    debug_file.write(f"\nИтого сделано 2 итерации\n")
    debug_file.write(f"Итоговое решение x = {x2.flatten().tolist()}\n")
    debug_file.close()
    
    print("\nОтладочная информация сохранена в файле py_debug.txt")
    print("Для сравнения с C++ версией, раскомментируйте строку с #define deb в msg_solver.cpp")
