# CE7453 Numerical Algorithms 期末高分超详细攻略（第六章：Eigenanalysis）

> **本章内容**：特征值与特征向量、幂迭代、QR算法、PageRank  
> **适用对象**：零基础/考前冲刺/快速查漏补缺  
> **内容特色**：详细原理、公式推导、算法流程、例题全解、英文关键词

---

## 1. 基本概念与背景

- **特征值与特征向量（Eigenvalue & Eigenvector）**：对于 $n \times n$ 矩阵 $A$，如果存在非零向量 $x$ 和标量 $\lambda$ 使 $Ax = \lambda x$，则 $\lambda$ 为特征值，$x$ 为对应特征向量。
- **实际应用**：PCA主成分分析、PageRank、振动分析、图像处理、机器学习等。

---

## 2. 幂迭代法（Power Iteration）

### 2.1 原理

- 用于求最大模特征值及其特征向量。
- 反复用 $A$ 作用于初始向量，逐步放大主特征分量。

### 2.2 算法流程

1. 选初始向量 $x_0$（非零）
2. 归一化 $x_0$
3. 迭代：$x_{k+1} = Ax_k$
4. 归一化 $x_{k+1}$
5. 若 $||x_{k+1} - x_k|| < \epsilon$，停止
6. 特征值近似为 $\lambda \approx \frac{x_k^T A x_k}{x_k^T x_k}$

### 2.3 例题

**例题1**：$A = \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}$，$x_0 = \begin{bmatrix} 1 \\ 0 \end{bmatrix}$，用幂迭代法求最大特征值

**解答**：
- $x_1 = A x_0 = [2,1]^T$，归一化 $[0.894, 0.447]^T$
- $x_2 = A x_1 = [2*0.894+1*0.447, 1*0.894+2*0.447] = [2.236, 1.788]$，归一化 $[0.781, 0.625]$
- 继续迭代，收敛到主特征向量 $[1,1]^T$，特征值约为 $3$

---

## 3. QR算法（QR Algorithm）

### 3.1 原理

- 用于求所有特征值（尤其适合对称矩阵）
- 反复分解 $A_k = Q_k R_k$，再 $A_{k+1} = R_k Q_k$

### 3.2 算法流程

1. 设 $A_0 = A$
2. 对 $A_k$ 做 QR分解：$A_k = Q_k R_k$
3. $A_{k+1} = R_k Q_k$
4. 重复，直到 $A_k$ 近似对角阵，对角线即为特征值

---

## 4. PageRank 算法（Google 搜索排序）

### 4.1 背景

- PageRank 本质是求网页链接矩阵的主特征向量
- 设 $P$ 为转移概率矩阵，PageRank 向量 $r$ 满足 $Pr = r$

### 4.2 算法流程

1. 构造转移矩阵 $P$
2. 选初始 $r_0$
3. 迭代 $r_{k+1} = P r_k$
4. 归一化 $r_{k+1}$
5. 收敛后 $r$ 即为网页排名

---

## 5. 常见考点与易错点

- 幂迭代初始向量不能为特征向量正交分量
- QR算法只适合小规模或对称矩阵
- 特征值/向量定义与实际意义
- PageRank 的矩阵构造与归一化
- 英文术语拼写与公式记忆

---

## 6. 英文关键词与表达

- eigenvalue, eigenvector, power iteration, Rayleigh quotient, QR algorithm, diagonalization, convergence, PageRank, principal component analysis (PCA), dominant eigenvalue

---

如需更详细例题推导或某一算法的代码实现，可随时补充！