# CE7453 Numerical Algorithms 期末高分超详细攻略（第三章：Bezier/B-spline/Interpolation）

> **本章内容**：Bezier 曲线、B-spline 曲线、插值方法  
> **适用对象**：零基础/考前冲刺/快速查漏补缺  
> **内容特色**：详细原理、公式推导、算法流程、例题全解、英文关键词

---

## 1. 基本概念与背景

- **Bezier 曲线**：由一组控制点定义的参数曲线，广泛用于计算机图形、字体、动画等。
- **B-spline 曲线**：Bezier 曲线的推广，支持更多控制点和更高阶的连续性，局部性更好。
- **插值（Interpolation）**：通过已知数据点，构造通过这些点的函数（曲线/多项式）。

---

## 2. Bezier 曲线

### 2.1 定义与公式

- **n阶 Bezier 曲线公式**：  
  $$
  P(t) = \sum_{i=0}^n B_{i,n}(t) P_i
  $$
  其中 $P_i$ 为控制点，$t \in [0,1]$，$B_{i,n}(t)$ 为 Bernstein 基函数：
  $$
  B_{i,n}(t) = \binom{n}{i} (1-t)^{n-i} t^i
  $$

### 2.2 主要性质

- **端点插值（Endpoint interpolation）**：$P(0)=P_0$，$P(1)=P_n$
- **共线性（Co-tangency）**：曲线在端点处与首末两段控制多边形共线
- **凸包性（Convex hull property）**：曲线始终在控制多边形的凸包内
- **仿射不变性（Affine invariance）**：对控制点做仿射变换，曲线同样变换

### 2.3 de Casteljau 算法（数值稳定的递归算法）

#### 算法流程
1. 设 $n+1$ 个控制点 $P_0, P_1, ..., P_n$，参数 $t \in [0,1]$
2. 递归计算：  
   $$
   P_i^{(0)} = P_i \\
   P_i^{(r)} = (1-t)P_i^{(r-1)} + tP_{i+1}^{(r-1)}, \quad r=1,2,...,n
   $$
3. 最终 $P_0^{(n)}$ 即为 $P(t)$

#### 步骤解释
- 每一层递归都在相邻点之间做线性插值（Lerp），逐步逼近曲线上的点。
- 该算法数值稳定，适合实际计算。

#### 例题
**已知控制点 $P_0=(0,0)$, $P_1=(1,2)$, $P_2=(3,3)$，求 $t=0.5$ 时的 Bezier 曲线点。**

**解答：**
- $P_0^{(0)}=(0,0)$, $P_1^{(0)}=(1,2)$, $P_2^{(0)}=(3,3)$
- $P_0^{(1)} = (1-0.5)(0,0) + 0.5(1,2) = (0.5,1)$  
  $P_1^{(1)} = (1-0.5)(1,2) + 0.5(3,3) = (2,2.5)$
- $P_0^{(2)} = (1-0.5)(0.5,1) + 0.5(2,2.5) = (1.25,1.75)$
- 答：$P(0.5) = (1.25, 1.75)$

**Python代码示例 (de Casteljau)**：
```python
import numpy as np

def de_casteljau(control_points, t):
    """使用 de Casteljau 算法计算 Bezier 曲线上一点"""
    points = np.array(control_points, dtype=float)
    n = len(points) - 1 # 曲线阶数
    
    # P_i^(0) = P_i
    current_points = points.copy()
    
    for r in range(1, n + 1): # 从 r=1 到 n
        new_points = []
        for i in range(n - r + 1): # P_i^(r) 的 i 从 0 到 n-r
            p_r_i = (1 - t) * current_points[i] + t * current_points[i+1]
            new_points.append(p_r_i)
        current_points = np.array(new_points)
        # print(f"r={r}, points={current_points}") # 可选：打印中间步骤
            
    return current_points[0] # P_0^(n)

# 例题数据
P_bezier = [[0, 0], [1, 2], [3, 3]]
t_bezier = 0.5

# 计算
point_on_curve = de_casteljau(P_bezier, t_bezier)
print(f"例题 (Bezier de Casteljau): P({t_bezier}) = {point_on_curve}")

# 验证 (可选)
# expected_point = np.array([1.25, 1.75])
# assert np.allclose(point_on_curve, expected_point)
```

---

## 3. B-spline 曲线

### 3.1 定义与公式

- **B-spline 曲线公式**：  
  $$
  r(u) = \sum_{i=0}^n N_{i,k}(u) P_i
  $$
  其中 $N_{i,k}(u)$ 为 $k$ 阶 B-spline 基函数，$P_i$ 为 de Boor 控制点，$u$ 为参数。

- **B-spline 基函数递归定义**：  
  $$
  N_{i,1}(u) = 
  \begin{cases}
    1, & u_i \leq u < u_{i+1} \\
    0, & \text{otherwise}
  \end{cases}
  $$
  $$
  N_{i,k}(u) = \frac{u-u_i}{u_{i+k-1}-u_i}N_{i,k-1}(u) + \frac{u_{i+k}-u}{u_{i+k}-u_{i+1}}N_{i+1,k-1}(u)
  $$

### 3.2 主要性质

- **局部性（Local support）**：每个控制点只影响相邻 $k$ 个区间
- **高阶连续性**：可实现 $C^{k-1}$ 连续
- **可调节点向量（Knot vector）**：灵活控制曲线形状

### 3.3 de Boor 算法

#### 算法流程
1. 给定节点向量 $U$，控制点 $P_0,...,P_n$，阶数 $k$，参数 $u$
2. 找到 $u$ 所在区间 $[u_j, u_{j+1})$
3. 递归线性插值，最终得到 $r(u)$

#### 例题
**已知三次 B-spline 曲线，控制点 $P_0=(0,0)$, $P_1=(1,2)$, $P_2=(3,3)$, $P_3=(4,0)$，节点向量 $[0,0,0,0,1,2,2,2,2]$，求 $u=1$ 处的曲线点。**

**解答：**
- 节点向量为 $[0,0,0,0,1,2,2,2,2]$，阶数 $k=4$（因为节点向量长度为 $n+k+1=4+4+1=9$）。
- 参数 $u=1$ 落在区间 $[u_3, u_4)=[0,1)$，但由于节点向量重复，实际有效区间需要考虑 $u=1$ 属于 $[u_4, u_5)=[1,2)$，因此 $j=4$。
- 影响 $u=1$ 的控制点为 $P_1, P_2, P_3$（因为 $j-k+1=4-4+1=1$ 到 $j=4$，即控制点索引 $1$ 到 $3$）。
- 计算基函数 $N_{i,4}(1)$：
  - 先计算 $N_{i,1}(1)$：
    - $N_{1,1}(1)=0$（因为 $u=1$ 不在 $[u_1,u_2)=[0,0)$）
    - $N_{2,1}(1)=0$（因为 $u=1$ 不在 $[u_2,u_3)=[0,0)$）
    - $N_{3,1}(1)=0$（因为 $u=1$ 不在 $[u_3,u_4)=[0,1)$）
    - $N_{4,1}(1)=1$（因为 $u=1$ 在 $[u_4,u_5)=[1,2)$）
    - $N_{5,1}(1)=0$（因为 $u=1$ 不在 $[u_5,u_6)=[2,2)$）
  - 递归计算 $N_{i,2}(1)$：
    - $N_{2,2}(1) = \frac{1-u_2}{u_3-u_2}N_{2,1}(1) + \frac{u_4-1}{u_4-u_3}N_{3,1}(1) = \frac{1-0}{0-0}*0 + \frac{1-1}{1-0}*0 = 0$
    - $N_{3,2}(1) = \frac{1-u_3}{u_4-u_3}N_{3,1}(1) + \frac{u_5-1}{u_5-u_4}N_{4,1}(1) = \frac{1-0}{1-0}*0 + \frac{2-1}{2-1}*1 = 1$
    - $N_{4,2}(1) = \frac{1-u_4}{u_5-u_4}N_{4,1}(1) + \frac{u_6-1}{u_6-u_5}N_{5,1}(1) = \frac{1-1}{2-1}*1 + \frac{2-1}{2-2}*0 = 0$
  - 递归计算 $N_{i,3}(1)$：
    - $N_{2,3}(1) = \frac{1-u_2}{u_4-u_2}N_{2,2}(1) + \frac{u_5-1}{u_5-u_3}N_{3,2}(1) = \frac{1-0}{1-0}*0 + \frac{2-1}{2-0}*1 = 0.5$
    - $N_{3,3}(1) = \frac{1-u_3}{u_5-u_3}N_{3,2}(1) + \frac{u_6-1}{u_6-u_4}N_{4,2}(1) = \frac{1-0}{2-0}*1 + \frac{2-1}{2-1}*0 = 0.5$
    - $N_{4,3}(1) = 0$（类似计算）
  - 递归计算 $N_{i,4}(1)$：
    - $N_{1,4}(1) = \frac{1-u_1}{u_4-u_1}N_{1,3}(1) + \frac{u_5-1}{u_5-u_2}N_{2,3}(1) = \frac{1-0}{1-0}*0 + \frac{2-1}{2-0}*0.5 = 0.25$
    - $N_{2,4}(1) = \frac{1-u_2}{u_5-u_2}N_{2,3}(1) + \frac{u_6-1}{u_6-u_3}N_{3,3}(1) = \frac{1-0}{2-0}*0.5 + \frac{2-1}{2-0}*0.5 = 0.5$
    - $N_{3,4}(1) = \frac{1-u_3}{u_6-u_3}N_{3,3}(1) + \frac{u_7-1}{u_7-u_4}N_{4,3}(1) = \frac{1-0}{2-0}*0.5 + \frac{2-1}{2-1}*0 = 0.25$
- 最终曲线点 $r(1) = N_{1,4}(1)P_1 + N_{2,4}(1)P_2 + N_{3,4}(1)P_3 = 0.25*(1,2) + 0.5*(3,3) + 0.25*(4,0) = (0.25*1 + 0.5*3 + 0.25*4, 0.25*2 + 0.5*3 + 0.25*0) = (2.75, 2.0)$
- 答：$r(1) = (2.75, 2.0)$

**Python代码示例 (手动计算基函数)**：
```python
import numpy as np

# 例题数据
control_points_bspline = np.array([[0, 0], [1, 2], [3, 3], [4, 0]])
knots_bspline = np.array([0,0,0,0,1,2,2,2,2])
u_bspline = 1

# 按照文档手动计算的基函数值
N14_at_1 = 0.25
N24_at_1 = 0.5
N34_at_1 = 0.25

# 影响 u=1 的控制点是 P1, P2, P3 (索引1, 2, 3)
point_on_bspline = N14_at_1 * control_points_bspline[1] + \
                   N24_at_1 * control_points_bspline[2] + \
                   N34_at_1 * control_points_bspline[3]
                   
print(f"例题 (B-spline 手动基函数): r({u_bspline}) = {point_on_bspline}")

# 验证 (可选)
# expected_point = np.array([2.75, 2.0])
# assert np.allclose(point_on_bspline, expected_point)

# 注意: 实际应用中会使用递归函数计算基函数或直接使用 scipy.interpolate
# from scipy.interpolate import BSpline
# k = 3 # B样条阶数 (degree), k=order-1. Order k=4 for cubic.
# tck = (knots_bspline, control_points_bspline.T, k)
# point_scipy = BSpline(*tck)(u_bspline)
# print(f"例题 (B-spline Scipy): r({u_bspline}) = {point_scipy}")
```

### 新增例题：计算给定控制点 $P_0=(0,0)$，$P_1=(1,1)$，$P_2=(2,0)$ 的二次Bezier曲线在 $t=0.5$ 处的点

**解答步骤**：
1. 使用de Casteljau算法：
   - 第一层：$P_0^{(1)} = 0.5 \times P_0 + 0.5 \times P_1 = (0.5, 0.5)$
   - $P_1^{(1)} = 0.5 \times P_1 + 0.5 \times P_2 = (1.5, 0.5)$
2. 第二层：$P_0^{(2)} = 0.5 \times P_0^{(1)} + 0.5 \times P_1^{(1)} = (1, 0.5)$
3. 结果：曲线点为 $(1, 0.5)$。

**Python代码示例 (使用上面的 de_casteljau 函数)**：
```python
# 新增例题数据
P_bezier_add = [[0, 0], [1, 1], [2, 0]]
t_bezier_add = 0.5

# 计算
point_on_curve_add = de_casteljau(P_bezier_add, t_bezier_add)
print(f"新增例题 (Bezier de Casteljau): P({t_bezier_add}) = {point_on_curve_add}")

# 验证 (可选)
# expected_point = np.array([1, 0.5])
# assert np.allclose(point_on_curve_add, expected_point)
```

---

## 4. 插值方法（Interpolation）

### 4.1 Lagrange 插值

- **公式**：  
  $$
  L(x) = \sum_{i=0}^n y_i \prod_{j\neq i} \frac{x-x_j}{x_i-x_j}
  $$
- **优缺点**：公式显式，适合手算，但高阶时数值不稳定。

#### 例题
**已知点 $(0,1)$, $(1,2)$, $(2,0)$，用 Lagrange 插值求 $x=1.5$ 处的函数值。**

**解答：**
- 基函数：
  - $l_0(x) = \frac{(x-1)(x-2)}{(0-1)(0-2)} = \frac{(x-1)(x-2)}{2}$
  - $l_1(x) = \frac{(x-0)(x-2)}{(1-0)(1-2)} = \frac{x(x-2)}{-1} = -x(x-2)$
  - $l_2(x) = \frac{(x-0)(x-1)}{(2-0)(2-1)} = \frac{x(x-1)}{2}$
- 插值多项式：$L(x) = 1*l_0(x) + 2*l_1(x) + 0*l_2(x) = \frac{(x-1)(x-2)}{2} - 2x(x-2)$
- 代入 $x=1.5$：
  - $l_0(1.5) = \frac{(1.5-1)(1.5-2)}{2} = \frac{0.5*(-0.5)}{2} = -0.125$
  - $l_1(1.5) = -1.5*(1.5-2) = -1.5*(-0.5) = 0.75$
  - $L(1.5) = 1*(-0.125) + 2*(0.75) = -0.125 + 1.5 = 1.375$
- 答：$L(1.5) = 1.375$

**Python代码示例 (Lagrange 插值)**：
```python
import numpy as np
from scipy import interpolate # 用于验证

def lagrange_basis(x_data, i, x):
    """计算第 i 个 Lagrange 基函数 l_i(x)"""
    n = len(x_data)
    li = 1.0
    for j in range(n):
        if i != j:
            li *= (x - x_data[j]) / (x_data[i] - x_data[j])
    return li

def lagrange_interpolation(x_data, y_data, x_interp):
    """使用 Lagrange 插值计算 x_interp 处的函数值"""
    n = len(x_data)
    result = 0.0
    for i in range(n):
        result += y_data[i] * lagrange_basis(x_data, i, x_interp)
    return result

# 例题数据
x_lagrange = np.array([0, 1, 2])
y_lagrange = np.array([1, 2, 0])
x_interp_lagrange = 1.5

# 计算
interp_value_lagrange = lagrange_interpolation(x_lagrange, y_lagrange, x_interp_lagrange)
print(f"例题 (Lagrange 插值): L({x_interp_lagrange}) = {interp_value_lagrange}")

# 验证 (可选)
# expected_value = 1.375
# assert np.isclose(interp_value_lagrange, expected_value)
# lagrange_poly = interpolate.lagrange(x_lagrange, y_lagrange)
# scipy_value = lagrange_poly(x_interp_lagrange)
# assert np.isclose(scipy_value, expected_value)
```

### 4.2 Newton 插值

- **公式**：  
  $$
  N(x) = a_0 + a_1(x-x_0) + a_2(x-x_0)(x-x_1) + \cdots
  $$
  其中 $a_i$ 为差商（divided differences）

#### 例题
**已知点 $(0,1)$, $(1,2)$, $(2,0)$，用 Newton 插值求 $x=1.5$ 处的函数值。**

**解答：**
- 差商表：
  - $f[0] = 1$
  - $f[1] = 2$
  - $f[2] = 0$
  - $f[0,1] = \frac{f[1]-f[0]}{1-0} = \frac{2-1}{1} = 1$
  - $f[1,2] = \frac{f[2]-f[1]}{2-1} = \frac{0-2}{1} = -2$
  - $f[0,1,2] = \frac{f[1,2]-f[0,1]}{2-0} = \frac{-2-1}{2} = -1.5$
- 插值多项式：$N(x) = 1 + 1*(x-0) + (-1.5)*(x-0)(x-1)$
- 代入 $x=1.5$：
  - $N(1.5) = 1 + 1*(1.5) + (-1.5)*(1.5)*(1.5-1) = 1 + 1.5 + (-1.5)*(1.5)*(0.5) = 1 + 1.5 - 1.125 = 1.375$
- 答：$N(1.5) = 1.375$

**Python代码示例 (Newton 插值)**：
```python
import numpy as np

def divided_differences(x_data, y_data):
    """计算差商表，返回 Newton 多项式系数"""
    n = len(x_data)
    coef = np.zeros([n, n])
    coef[:, 0] = y_data # 第一列是 y 值
    
    for j in range(1, n): # 列
        for i in range(n - j): # 行
            coef[i, j] = (coef[i + 1, j - 1] - coef[i, j - 1]) / (x_data[i + j] - x_data[i])
            
    return coef[0, :] # 返回第一行作为系数 a0, a1, a2, ...

def newton_polynomial(coefficients, x_data, x):
    """根据差商系数计算 Newton 插值多项式在 x 处的值"""
    n = len(coefficients)
    result = coefficients[0]
    term = 1.0
    
    for i in range(1, n):
        term *= (x - x_data[i - 1]) # 计算 (x-x0), (x-x0)(x-x1), ...
        result += coefficients[i] * term
        
    return result

# 例题数据
x_newton = np.array([0, 1, 2])
y_newton = np.array([1, 2, 0])
x_interp_newton = 1.5

# 计算
coeffs_newton = divided_differences(x_newton, y_newton)
interp_value_newton = newton_polynomial(coeffs_newton, x_newton, x_interp_newton)
print(f"例题 (Newton 插值) 差商系数: {coeffs_newton}")
print(f"例题 (Newton 插值): N({x_interp_newton}) = {interp_value_newton}")

# 验证 (可选)
# expected_value = 1.375
# assert np.isclose(interp_value_newton, expected_value)
```

### 4.3 样条插值（Spline Interpolation）

- **三次样条**：分段三次多项式，保证 $C^2$ 连续，常用于平滑曲线拟合。

#### 例题
**已知点 $(0,1)$, $(1,2)$, $(2,0)$，构造自然三次样条插值（边界条件为二阶导数为0）。**

**解答：**
- 设区间 $[0,1]$ 和 $[1,2]$ 上的三次多项式为：
  - $S_1(x) = a_1 + b_1(x-0) + c_1(x-0)^2 + d_1(x-0)^3$
  - $S_2(x) = a_2 + b_2(x-1) + c_2(x-1)^2 + d_2(x-1)^3$
- 条件：
  1. 插值条件：$S_1(0)=1$, $S_1(1)=2$, $S_2(1)=2$, $S_2(2)=0$
  2. 一阶导数连续：$S_1'(1)=S_2'(1)$
  3. 二阶导数连续：$S_1''(1)=S_2''(1)$
  4. 自然边界条件：$S_1''(0)=0$, $S_2''(2)=0$
- 解方程组（详细计算略）：
  - $S_1(x) = 1 + 1.5x - 0.5x^3$
  - $S_2(x) = 2 - 0.5(x-1) - 1.5(x-1)^2 + 0.5(x-1)^3$
- 答：自然三次样条插值多项式如上。

**Python代码示例 (使用 Scipy)**：
```python
import numpy as np
from scipy.interpolate import CubicSpline

# 例题数据
x_spline = np.array([0, 1, 2])
y_spline = np.array([1, 2, 0])

# 计算自然三次样条 (bc_type='natural')
cs = CubicSpline(x_spline, y_spline, bc_type='natural')

# 打印样条系数 (可选)
# print("样条系数 (分段):")
# for i in range(len(cs.c[0, :])):
#     print(f" 区间 {i}: {cs.c[:, i]}")

# 在 x=1.5 处插值 (虽然例题只要求构造)
x_interp_spline = 1.5
interp_value_spline = cs(x_interp_spline)
print(f"例题 (自然三次样条插值): S({x_interp_spline}) = {interp_value_spline}")

# 验证 x=0, 1, 2 处的插值结果
print(f"验证: S(0)={cs(0)}, S(1)={cs(1)}, S(2)={cs(2)}")
assert np.isclose(cs(0), 1)
assert np.isclose(cs(1), 2)
assert np.isclose(cs(2), 0)

# 验证边界条件 S''(0)=0, S''(2)=0
print(f"验证: S''(0)={cs(0, 2)}, S''(2)={cs(2, 2)}")
assert np.isclose(cs(0, 2), 0)
assert np.isclose(cs(2, 2), 0)
```

---

## 5. 常见考点与易错点

- Bezier 曲线的端点、切线、凸包等性质
- de Casteljau 算法步骤与递归思想
- B-spline 节点向量与控制点影响范围
- 插值多项式的数值稳定性与误差
- 英文术语拼写与公式记忆

---

## 6. 英文关键词与表达

- Bezier curve, control point, Bernstein polynomial, de Casteljau algorithm, subdivision, convex hull, affine invariance
- B-spline, knot vector, de Boor point, local support, continuity, degree, order
- interpolation, Lagrange interpolation, Newton divided difference, spline, cubic spline

---

如需更详细例题推导或某一算法的代码实现，可随时补充！