import numpy as np
from scipy.optimize import fsolve
import math
import matplotlib.pyplot as plt  # 导入 matplotlib 进行绘图
import matplotlib.patches as patches  # 导入 patches 以绘制矩形

# 定义常量
b = 286

# 定义需要解的方程
def equation(y, x, b):
    return 2 * math.cos(y) * (x * (x + 27.5 / np.pi * y)) - (x * x + (x + 27.5 / np.pi * y) * (x + 27.5 / np.pi * y) - b * b)

# 初始化 x_array
x_array = np.array([300])  # 这是一个包含初始值的 NumPy 一维数组

# 存储所有笛卡尔坐标
all_x_cartesian = []
all_y_cartesian = []

# 打开文件用于写入
with open('第三问极坐标.txt', 'w') as file:
    # 进行 1 次循环
    for iteration in range(1):
        # 存储解 y 的列表
        y_solutions = []

        # 对当前的 x_array 中的每个 x 值求解 y
        for x in x_array:
            # 使用 fsolve 来求解 y 的值，初始猜测 y = 1.0
            y_solution = fsolve(equation, 1.0, args=(x, b))
            y_solutions.append(y_solution[0])

        # 计算 z_array (极角数组)
        z_array = x_array * np.pi / 27.5

        # 将极坐标转换为笛卡尔坐标
        x_cartesian = x_array * np.cos(z_array)  # 极径 x_array 作为 r
        y_cartesian = x_array * np.sin(z_array)  # 极角 z_array 作为 θ

        # 存储到所有坐标列表
        all_x_cartesian.extend(x_cartesian)
        all_y_cartesian.extend(y_cartesian)

        # 输出当前循环的结果并写入文件
        file.write(f"第 {iteration + 1} 次循环结果:\n")
        for i in range(len(x_array)):
            file.write(f"极坐标: r = {x_array[i]:.6f}, θ = {z_array[i]:.6f} 弧度 -> 笛卡尔坐标: x = {x_cartesian[i]:.6f}, y = {y_cartesian[i]:.6f}\n")

        # 更新 x_array
        x_array = np.array(y_solutions) * 27.5 / np.pi + x_array

    # 定义常量
    b = 165

    # 使用上个循环的 x_array
    for iteration in range(15):
        # 存储解 y 的列表
        y_solutions = []

        # 对当前的 x_array 中的每个 x 值求解 y
        for x in x_array:
            y_solution = fsolve(equation, 1.0, args=(x, b))
            y_solutions.append(y_solution[0])

        # 计算 z_array (极角数组)
        z_array = x_array * np.pi / 27.5

        # 将极坐标转换为笛卡尔坐标
        x_cartesian = x_array * np.cos(z_array)  # 极径 x_array 作为 r
        y_cartesian = x_array * np.sin(z_array)  # 极角 z_array 作为 θ

        # 存储到所有坐标列表
        all_x_cartesian.extend(x_cartesian)
        all_y_cartesian.extend(y_cartesian)

        # 输出当前循环的结果并写入文件
        file.write(f"\n第 {iteration + 1} 次循环结果:\n")
        for i in range(len(x_array)):
            file.write(f"极坐标: r = {x_array[i]:.6f}, θ = {z_array[i]:.6f} 弧度 -> 笛卡尔坐标: x = {x_cartesian[i]:.6f}, y = {y_cartesian[i]:.6f}\n")

        # 更新 x_array
        x_array = np.array(y_solutions) * 27.5 / np.pi + x_array

print("结果已保存到 'output_data3.txt'")

# 可视化笛卡尔坐标并将点连起来
plt.figure(figsize=(8, 8))
plt.plot(all_x_cartesian, all_y_cartesian, marker='o', linestyle='-')  # 使用 plot 连线

plt.title("极坐标转换为笛卡尔坐标的可视化 (点与线连接) + 中点矩形")
plt.xlabel("X 坐标")
plt.ylabel("Y 坐标")
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()





