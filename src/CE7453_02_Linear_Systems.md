# CE7453 Numerical Algorithms 期末高分超详细攻略（第二章：Linear Systems）

> **本章内容**：线性方程组（Linear Systems）  
> **适用对象**：零基础/考前冲刺/快速查漏补缺  
> **内容特色**：详细原理、公式推导、算法流程、例题全解、英文关键词

---

## 2.1 基本概念与背景

- **什么是线性方程组？**  
  一组形如 $Ax=b$ 的方程，$A$ 是矩阵，$x$ 是未知向量，$b$ 是常数向量。
- **实际应用**：  
  工程、物理、数据拟合、网络流等大量问题都可归结为解线性方程组。

---

## 2.2 主要算法与原理

### 2.2.1 高斯消元法（Gaussian Elimination）

- **原理**：  
  通过行变换，把矩阵化为上三角，再回代求解。
- **算法流程**：
  1. 对每一列，从当前行到最后一行找最大主元（pivoting，防止数值不稳定）。
  2. 用当前行消去下面各行的该列元素。
  3. 得到上三角矩阵后，从最后一行开始回代求解未知数。
- **公式**：  
  $a_{ij}^{(k+1)} = a_{ij}^{(k)} - m_{ik}a_{kj}^{(k)}$，$m_{ik}=a_{ik}^{(k)}/a_{kk}^{(k)}$

- **英文关键词**：Gaussian elimination, row operation, pivoting, upper triangular, back substitution

---

### 2.2.2 LU分解（LU Decomposition）

- **原理**：  
  把 $A$ 分解为 $LU$，$L$下三角，$U$上三角。  
  先解 $Ly=b$，再解 $Ux=y$。
- **优点**：  
  多次解 $Ax=b$ 时只需分解一次 $A$，高效。

- **英文关键词**：LU decomposition, lower triangular, upper triangular, forward substitution, backward substitution

---

### 2.2.3 Jacobi与Gauss-Seidel迭代法

- **Jacobi法公式**：  
  $x_i^{(k+1)} = \frac{1}{a_{ii}}(b_i - \sum_{j\neq i} a_{ij}x_j^{(k)})$
- **Gauss-Seidel法公式**：  
  $x_i^{(k+1)} = \frac{1}{a_{ii}}(b_i - \sum_{j<i} a_{ij}x_j^{(k+1)} - \sum_{j>i} a_{ij}x_j^{(k)})$
- **收敛条件**：  
  系数矩阵严格对角占优（diagonally dominant）时保证收敛。

- **英文关键词**：Jacobi method, Gauss-Seidel method, iteration, convergence, diagonal dominance

---

### 2.2.4 梯度下降法（Gradient Descent）

- **原理**：  
  对称正定矩阵 $A$，$Ax=b$ 可转化为最小化 $f(x)=\frac{1}{2}x^TAx-b^Tx$。
- **迭代公式**：  
  $x_{k+1} = x_k + \alpha_k r_k$，$r_k = b - Ax_k$，$\alpha_k = \frac{r_k^Tr_k}{r_k^TAr_k}$

- **英文关键词**：gradient descent, positive definite, minimization, residual

---

## 2.3 典型例题与详细解答

### 例题1：用高斯消元法解
$$
\begin{cases}
2x + 3y = 8 \\
5x + 4y = 13
\end{cases}
$$

**解答步骤**：
1. 用第一行消去第二行的 $x$：
   - $m = 5/2 = 2.5$
   - 第二行 $-2.5 \times$ 第一行：$[5,4,13] - 2.5*[2,3,8] = [0, -3.5, -7]$
2. 回代：$-3.5y = -7 \Rightarrow y=2$
3. 代入第一行：$2x+3*2=8 \Rightarrow x=1$

---

### 例题2：用Jacobi法迭代一次
$$
\begin{cases}
4x + y = 9 \\
x + 3y = 7
\end{cases}
$$
初始 $x^{(0)}=0, y^{(0)}=0$

**解答步骤**：
- $x^{(1)} = (9 - 0)/4 = 2.25$
- $y^{(1)} = (7 - 0)/3 = 2.333...$

---

## 2.4 常见考点与易错点

- 主元选取不当导致数值不稳定
- 迭代法收敛条件不满足
- 回代顺序错误
- 复杂度分析（高斯消元 $O(n^3)$）