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

**解答思路**：
- 找到 $u=1$ 所在区间
- 按 de Boor 算法递归插值
- 详细步骤略（实际考试可按讲义模板写出每步）

---

## 4. 插值方法（Interpolation）

### 4.1 Lagrange 插值

- **公式**：  
  $$
  L(x) = \sum_{i=0}^n y_i \prod_{j\neq i} \frac{x-x_j}{x_i-x_j}
  $$
- **优缺点**：公式显式，适合手算，但高阶时数值不稳定。

### 4.2 Newton 插值

- **公式**：  
  $$
  N(x) = a_0 + a_1(x-x_0) + a_2(x-x_0)(x-x_1) + \cdots
  $$
  其中 $a_i$ 为差商（divided differences）

### 4.3 样条插值（Spline Interpolation）

- **三次样条**：分段三次多项式，保证 $C^2$ 连续，常用于平滑曲线拟合。

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