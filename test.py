import numpy as np
from scipy.optimize import fsolve
import math
import matplotlib.pyplot as plt  # 导入 matplotlib 进行绘图

# 定义常量
b = 286

# 定义需要解的方程
def equation(y, x, b):
    return 2 * math.cos(y) * (x * (x + 21.6 / np.pi * y)) - (x * x + (x + 21.6 / np.pi * y) * (x + 21.6 / np.pi * y) - b * b)

# 初始化 x_array
x_array = np.array([450])  # 这是一个包含初始值的 NumPy 一维数组
print(x_array)
# 存储所有笛卡尔坐标
all_x_cartesian = []
all_y_cartesian = []
distances_left = []
distances_right = [] 
min_distances = []
distance_left = []
distance_right = []
min_distance = []
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
        z_array = x_array * np.pi / 21.6

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
        x_array = np.array(y_solutions) * 21.6 / np.pi + x_array

    # 定义常量
    b = 165

    # 使用上个循环的 x_array
    for iteration in range(50):
        # 存储解 y 的列表
        y_solutions = []

        # 对当前的 x_array 中的每个 x 值求解 y
        for x in x_array:
            y_solution = fsolve(equation, 1.0, args=(x, b))
            y_solutions.append(y_solution[0])

        # 计算 z_array (极角数组)
        z_array = x_array * np.pi / 21.6

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
        x_array = np.array(y_solutions) * 21.6 / np.pi + x_array

# 自动使用前两个笛卡尔坐标点
x1, y1 = all_x_cartesian[0], all_y_cartesian[0]
x2, y2 = all_x_cartesian[1], all_y_cartesian[1]

# 计算两个点之间的中点
def midpoint(x1, y1, x2, y2):
    return ( (x1 + x2) / 2, (y1 + y2) / 2 )

# 计算两点连线的方向向量，并生成单位方向向量
def direction_vector(x1, y1, x2, y2):
    dx = x2 - x1
    dy = y2 - y1
    length = np.sqrt(dx**2 + dy**2)
    return (dx / length, dy / length)

# 计算长方形的四个顶点
def rectangle_vertices(x1, y1, x2, y2, length=341, width=30):
    # 计算中点
    mx, my = midpoint(x1, y1, x2, y2)

    # 计算方向向量（单位向量）
    dx, dy = direction_vector(x1, y1, x2, y2)

    # 长方形的长边的一半
    half_length = length / 2
    # 长方形的宽的一半
    half_width = width / 2

    # 旋转90度方向向量，得到垂直方向的单位向量
    perp_dx = -dy
    perp_dy = dx

    # 计算两个外侧的顶点

    bottom_left_x = mx - half_length * dx - half_width * perp_dx
    bottom_left_y = my - half_length * dy - half_width * perp_dy

    bottom_right_x = mx + half_length * dx - half_width * perp_dx
    bottom_right_y = my + half_length * dy - half_width * perp_dy

    # 返回四个顶点的坐标
    return [(bottom_left_x, bottom_left_y), (bottom_right_x, bottom_right_y)]

# 计算长方形的顶点
vertices = rectangle_vertices(x1, y1, x2, y2)

# 输出两个顶点的坐标
print("龙头的两个外侧顶点的坐标是：")
for i, (x, y) in enumerate(vertices, start=1):
    print(f"顶点 {i}: ({x:.2f}, {y:.2f})")

# 计算点到直线的距离
def point_to_line_distance(x0, y0, x1, y1, x2, y2):
    # 点 (x0, y0) 到直线 (x1, y1)-(x2, y2) 的距离公式
    distance = abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
    return distance

# 计算 bottom_left 和 bottom_right 到从第三个点往后的所有线段的距离
for i in range(2, len(all_x_cartesian) - 1):
    line_x1, line_y1 = all_x_cartesian[i], all_y_cartesian[i]
    line_x2, line_y2 = all_x_cartesian[i + 1], all_y_cartesian[i + 1]

    # 计算 bottom_left 到当前线段的距离
    distance_left = point_to_line_distance(vertices[0][0], vertices[0][1], line_x1, line_y1, line_x2, line_y2)
    distances_left.append(distance_left)
    # 计算 bottom_right 到当前线段的距离
    distance_right = point_to_line_distance(vertices[1][0], vertices[1][1], line_x1, line_y1, line_x2, line_y2)
    distances_right.append(distance_right)
     # 计算并输出最小距离
min_distance_left = min(distances_left)
min_distance_right = min(distances_right)
min_distance = min(min_distance_left, min_distance_right)
print(f"最小距离: {min_distance:.2f}")