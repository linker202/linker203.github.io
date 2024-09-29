import math

def calculate_distance(x1, y1, x2, y2):
    # 使用欧几里得距离公式计算两点之间的距离
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return distance

# 示例：输入两点的坐标
x1, y1 = 8.73 , 0.01 # 第一点坐标
x2, y2 = 8.62 , 2.86 # 第二点坐标

# 调用函数计算距离
distance = calculate_distance(x1, y1, x2, y2)
print(f"点 ({x1}, {y1}) 和点 ({x2}, {y2}) 之间的距离是: {distance:.2f}")