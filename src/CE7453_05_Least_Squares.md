# CE7453 Numerical Algorithms 期末高分超详细攻略（第五章：Least Squares）

> **本章内容**：最小二乘法、正规方程、QR分解、非线性最小二乘、GPS应用  
> **适用对象**：零基础/考前冲刺/快速查漏补缺  
> **内容特色**：详细原理、公式推导、算法流程、例题全解、英文关键词

---

## 1. 基本概念与背景

- **最小二乘法（Least Squares）**：用于拟合数据、解超定方程组（方程数多于未知数），使误差平方和最小。
- **实际应用**：数据拟合、回归分析、信号处理、GPS定位、机器学习等。

---

## 2. 线性最小二乘法

### 2.1 问题描述

- 给定 $m$ 个方程 $Ax \approx b$，$A$ 为 $m \times n$ 矩阵（$m>n$），通常无精确解。
- 目标：找到 $x$ 使 $||Ax-b||^2$ 最小。

### 2.2 正规方程（Normal Equation）

- **推导**：
  $$
  \min_x ||Ax-b||^2 = (Ax-b)^T(Ax-b)
  $$
  对 $x$ 求导，令导数为0，得
  $$
  A^TAx = A^Tb
  $$
- **解法**：解 $n$ 阶方程组 $A^TAx = A^Tb$

### 2.3 QR分解法（QR Decomposition）

- **原理**：将 $A$ 分解为 $A=QR$，$Q$ 为正交矩阵，$R$ 为上三角矩阵。
- **步骤**：
  1. $A=QR$
  2. $Rx=Q^Tb$，回代求解 $x$
- **优点**：数值稳定性高，适合大规模问题

### 2.4 典型例题

**例题1**：用最小二乘法拟合直线 $y=ax+b$，已知数据点 $(1,2), (2,3), (3,5)$

**解答**：
- 写成矩阵形式：
  $$
  \begin{bmatrix}
  1 & 1 \\
  2 & 1 \\
  3 & 1
  \end{bmatrix}
  \begin{bmatrix}
  a \\ b
  \end{bmatrix}
  =
  \begin{bmatrix}
  2 \\ 3 \\ 5
  \end{bmatrix}
  $$
- $A^TA = \begin{bmatrix} 14 & 6 \\ 6 & 3 \end{bmatrix}$，$A^Tb = \begin{bmatrix} 23 \\ 10 \end{bmatrix}$
- 解 $A^TAx = A^Tb$，得 $a=1.5, b=0.333...$

---

## 3. 非线性最小二乘（Nonlinear Least Squares）

### 3.1 问题描述

- 拟合模型 $y=f(x,\theta)$，$\theta$ 为参数，$f$ 非线性
- 目标：最小化 $S(\theta) = \sum_{i=1}^m (y_i - f(x_i, \theta))^2$

### 3.2 Gauss-Newton 方法

- **思想**：用泰勒展开线性化，迭代求解
- **步骤**：
  1. 初始猜测 $\theta_0$
  2. 计算残差 $r_i = y_i - f(x_i, \theta)$
  3. 计算雅可比矩阵 $J_{ij} = \frac{\partial f(x_i, \theta)}{\partial \theta_j}$
  4. 解正规方程 $J^TJ \Delta\theta = J^Tr$
  5. 更新 $\theta \leftarrow \theta + \Delta\theta$
  6. 迭代至收敛

---

## 4. GPS定位中的最小二乘

- **原理**：通过测量与多颗卫星的距离，建立超定方程组，利用最小二乘法估算接收机位置。
- **模型**：$||x - s_i|| = d_i$，$x$为未知位置，$s_i$为卫星位置，$d_i$为测距
- **步骤**：线性化后用Gauss-Newton迭代

---

## 5. 常见考点与易错点

- 正规方程推导与矩阵乘法
- QR分解步骤与正交化方法（Gram-Schmidt）
- 非线性最小二乘的迭代过程
- 残差（residual）、条件数（condition number）理解
- 英文术语拼写与公式记忆

---

## 6. 英文关键词与表达

- least squares, overdetermined system, normal equation, QR decomposition, orthogonalization, Gram-Schmidt, residual, regression, nonlinear least squares, Gauss-Newton, GPS positioning, condition number

---

如需更详细例题推导或某一算法的代码实现，可随时补充！