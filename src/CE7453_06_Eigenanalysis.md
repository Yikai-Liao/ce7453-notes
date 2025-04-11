# CE7453 Numerical Algorithms 期末高分超详细攻略（第六章：Eigenanalysis）

> **本章内容**：特征值与特征向量、幂迭代、QR算法、PageRank  
> **适用对象**：零基础/考前冲刺/快速查漏补缺  
> **内容特色**：详细原理、公式推导、算法流程、例题全解、英文关键词

---

## 1. 基本概念与背景

- **特征值与特征向量（Eigenvalue & Eigenvector）**：对于 $n \times n$ 矩阵 $A$，如果存在非零向量 $x$ 和标量 $\lambda$ 使 $Ax = \lambda x$，则 $\lambda$ 为特征值，$x$ 为对应特征向量。
- **实际应用**：PCA主成分分析、PageRank、振动分析、图像处理、机器学习等。
- **数学意义**：特征值表示矩阵在某些方向上的缩放因子，特征向量表示这些方向。特征分解可以帮助理解矩阵的性质和行为。

**特征值分解（Eigenvalue Decomposition）**：
如果矩阵 $A$ 可对角化，则可以表示为 $A = Q \Lambda Q^{-1}$，其中 $Q$ 是特征向量组成的矩阵，$\Lambda$ 是特征值组成对角矩阵。这种分解在求矩阵幂、解微分方程等场景中非常有用。

---

## 2. 幂迭代法（Power Iteration）

### 2.1 原理

- 用于求最大模特征值（dominant eigenvalue）及其对应的特征向量。
- 基本思想是通过反复用矩阵 $A$ 作用于初始向量，逐步放大主特征分量（与最大特征值对应的分量），最终收敛到主特征向量。
- 数学依据：假设矩阵 $A$ 有特征值 $\lambda_1, \lambda_2, ..., \lambda_n$，且 $|\lambda_1| > |\lambda_2| \geq ... \geq |\lambda_n|$，初始向量 $x_0$ 可以表示为特征向量的线性组合 $x_0 = c_1 v_1 + c_2 v_2 + ... + c_n v_n$，则 $A^k x_0 \approx c_1 \lambda_1^k v_1$（当 $k$ 很大时），即收敛到主特征向量方向。

### 2.2 算法流程

1. 选择一个非零初始向量 $x_0$（通常随机选择或全为1）。
2. 归一化 $x_0$，以避免数值过大或过小。
3. 迭代：计算 $x_{k+1} = A x_k$。
4. 归一化 $x_{k+1}$，得到单位向量。
5. 检查收敛性：若 $||x_{k+1} - x_k|| < \epsilon$（或达到最大迭代次数），则停止迭代。
6. 计算特征值近似：$\lambda \approx \frac{x_k^T A x_k}{x_k^T x_k}$（Rayleigh quotient，瑞利商）。

**注意事项**：
- 初始向量不能与主特征向量正交，否则无法收敛到最大特征值。
- 收敛速度取决于特征值之间的比值 $|\lambda_2 / \lambda_1|$，比值越小，收敛越快。

### 2.3 例题

**例题1**：给定矩阵 $A = \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}$，初始向量 $x_0 = \begin{bmatrix} 1 \\ 0 \end{bmatrix}$，用幂迭代法求最大特征值及其特征向量。

**详细解答**：
1. **初始化**：$x_0 = \begin{bmatrix} 1 \\ 0 \end{bmatrix}$，归一化后仍为 $\begin{bmatrix} 1 \\ 0 \end{bmatrix}$（因为 $||x_0||_2 = 1$）。
2. **第一次迭代**：
   - 计算 $x_1 = A x_0 = \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix} \begin{bmatrix} 1 \\ 0 \end{bmatrix} = \begin{bmatrix} 2 \\ 1 \end{bmatrix}$。
   - 归一化：$||x_1||_2 = \sqrt{2^2 + 1^2} = \sqrt{5} \approx 2.236$，所以 $x_1 = \begin{bmatrix} 2/2.236 \\ 1/2.236 \end{bmatrix} \approx \begin{bmatrix} 0.894 \\ 0.447 \end{bmatrix}$。
3. **第二次迭代**：
   - 计算 $x_2 = A x_1 = \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix} \begin{bmatrix} 0.894 \\ 0.447 \end{bmatrix} = \begin{bmatrix} 2*0.894 + 1*0.447 \\ 1*0.894 + 2*0.447 \end{bmatrix} = \begin{bmatrix} 2.235 \\ 1.788 \end{bmatrix}$。
   - 归一化：$||x_2||_2 = \sqrt{2.235^2 + 1.788^2} \approx 2.86$，所以 $x_2 \approx \begin{bmatrix} 0.781 \\ 0.625 \end{bmatrix}$。
4. **第三次迭代**：
   - 计算 $x_3 = A x_2 \approx \begin{bmatrix} 2*0.781 + 1*0.625 \\ 1*0.781 + 2*0.625 \end{bmatrix} = \begin{bmatrix} 2.187 \\ 2.031 \end{bmatrix}$。
   - 归一化：$||x_3||_2 \approx 2.98$，所以 $x_3 \approx \begin{bmatrix} 0.734 \\ 0.681 \end{bmatrix}$。
5. **第四次迭代**：
   - 计算 $x_4 = A x_3 \approx \begin{bmatrix} 2*0.734 + 1*0.681 \\ 1*0.734 + 2*0.681 \end{bmatrix} = \begin{bmatrix} 2.149 \\ 2.096 \end{bmatrix}$。
   - 归一化：$||x_4||_2 \approx 3.00$，所以 $x_4 \approx \begin{bmatrix} 0.716 \\ 0.698 \end{bmatrix}$。
6. **第五次迭代**：
   - 计算 $x_5 = A x_4 \approx \begin{bmatrix} 2*0.716 + 1*0.698 \\ 1*0.716 + 2*0.698 \end{bmatrix} = \begin{bmatrix} 2.130 \\ 2.112 \end{bmatrix}$。
   - 归一化：$||x_5||_2 \approx 3.00$，所以 $x_5 \approx \begin{bmatrix} 0.710 \\ 0.704 \end{bmatrix}$。
7. **收敛检查**：继续迭代，向量逐渐接近 $\begin{bmatrix} 0.707 \\ 0.707 \end{bmatrix}$（即 $\begin{bmatrix} 1 \\ 1 \end{bmatrix}$ 归一化后），说明收敛到主特征向量。
8. **计算特征值**：使用瑞利商，$\lambda \approx \frac{x_5^T A x_5}{x_5^T x_5}$。
   - $A x_5 \approx \begin{bmatrix} 2*0.710 + 1*0.704 \\ 1*0.710 + 2*0.704 \end{bmatrix} = \begin{bmatrix} 2.124 \\ 2.118 \end{bmatrix}$。
   - $x_5^T A x_5 \approx 0.710*2.124 + 0.704*2.118 \approx 3.000$。
   - $x_5^T x_5 = 0.710^2 + 0.704^2 \approx 1$。
   - 因此 $\lambda \approx 3$。

**答案**：最大特征值约为 $\lambda = 3$，对应特征向量为 $\begin{bmatrix} 1 \\ 1 \end{bmatrix}$（未归一化形式）。

**验证**：直接计算特征值，$A$ 的特征方程为 $\det(A - \lambda I) = (2-\lambda)^2 - 1 = 0$，解得 $\lambda = 3$ 或 $\lambda = 1$，与幂迭代结果一致。

---

## 3. QR算法（QR Algorithm）

### 3.1 原理

- QR算法是一种求解矩阵所有特征值的迭代方法，特别适合对称矩阵。
- 基本思想：通过反复对矩阵进行QR分解（将矩阵分解为正交矩阵 $Q$ 和上三角矩阵 $R$ 的乘积），并重新组合为 $A_{k+1} = R_k Q_k$，最终使矩阵收敛到对角形式（或近似对角形式），对角线元素即为特征值。
- 数学依据：QR算法本质上是幂法的扩展，同时对所有特征值进行迭代，收敛后矩阵变为上三角或对角矩阵。

### 3.2 算法流程

1. 初始化：令 $A_0 = A$。
2. 对当前矩阵 $A_k$ 进行QR分解：$A_k = Q_k R_k$（其中 $Q_k$ 是正交矩阵，$R_k$ 是上三角矩阵）。
3. 计算新的矩阵：$A_{k+1} = R_k Q_k$。
4. 重复步骤2和3，直到 $A_k$ 近似为对角矩阵（或达到最大迭代次数），对角线上的元素即为特征值的近似值。

**注意事项**：
- 对于对称矩阵，QR算法会收敛到对角矩阵；对于非对称矩阵，收敛到上三角矩阵（Schur形式）。
- 实际应用中常结合Hessenberg化（将矩阵化为上Hessenberg形式）和位移策略（shift）来加速收敛。

### 3.3 例题

**例题2**：使用QR算法计算矩阵 $A = \begin{bmatrix} 5 & 4 \\ 4 & 5 \end{bmatrix}$ 的特征值。

**详细解答**：
1. **初始化**：$A_0 = \begin{bmatrix} 5 & 4 \\ 4 & 5 \end{bmatrix}$。
2. **第一次迭代**：
   - 对 $A_0$ 进行QR分解：
     - 使用Gram-Schmidt正交化或其他方法，得到 $Q_0 = \begin{bmatrix} 0.707 & -0.707 \\ 0.707 & 0.707 \end{bmatrix}$，$R_0 = \begin{bmatrix} 7.071 & 6.364 \\ 0 & 1.414 \end{bmatrix}$（近似值）。
   - 计算 $A_1 = R_0 Q_0 = \begin{bmatrix} 7.071 & 6.364 \\ 0 & 1.414 \end{bmatrix} \begin{bmatrix} 0.707 & -0.707 \\ 0.707 & 0.707 \end{bmatrix} \approx \begin{bmatrix} 9.5 & -0.5 \\ 1.0 & 0.5 \end{bmatrix}$。
3. **第二次迭代**：
   - 对 $A_1$ 进行QR分解，得到新的 $Q_1$ 和 $R_1$。
   - 计算 $A_2 = R_1 Q_1$，结果更接近对角矩阵。
4. **继续迭代**：经过多次迭代，矩阵逐渐收敛到 $\begin{bmatrix} 9 & 0 \\ 0 & 1 \end{bmatrix}$，即特征值为 $9$ 和 $1$。

**答案**：特征值为 $\lambda_1 = 9$，$\lambda_2 = 1$。

**验证**：直接计算特征方程 $\det(A - \lambda I) = (5-\lambda)^2 - 16 = 0$，解得 $\lambda = 9$ 或 $\lambda = 1$，与QR算法结果一致。

---

## 4. PageRank 算法（Google 搜索排序）

### 4.1 背景

- PageRank 是Google搜索引擎的核心算法之一，用于评估网页的重要性。
- 基本思想：网页的重要性取决于链接到该网页的其他网页数量和重要性，本质上是求网页链接矩阵的主特征向量。
- 数学模型：设 $P$ 为转移概率矩阵（网页之间的链接概率），PageRank向量 $r$ 满足 $P r = r$，即 $r$ 是矩阵 $P$ 对应特征值 $\lambda=1$ 的特征向量。

### 4.2 算法流程

1. 构造转移矩阵 $P$：$P_{ij}$ 表示从网页 $j$ 链接到网页 $i$ 的概率（若网页 $j$ 有 $k$ 个外链，则每个外链概率为 $1/k$）。
2. 选择初始向量 $r_0$（通常为均匀分布，即每个网页初始重要性相等）。
3. 迭代：$r_{k+1} = P r_k$。
4. 归一化 $r_{k+1}$，确保向量元素之和为1。
5. 收敛后，$r$ 的各个分量即为各网页的PageRank值（排名分数）。

**阻尼因子（Damping Factor）**：
- 实际中引入阻尼因子 $d$（通常为0.85），以模拟用户随机跳转行为，修正公式为 $r = (1-d)/n + d P r$，其中 $n$ 为网页总数。
- 这保证了矩阵的收敛性，避免某些网页没有外链导致的问题。

### 4.3 例题

**例题3**：假设有3个网页，链接关系如下：网页1链接到2和3，网页2链接到3，网页3链接到1。计算各网页的PageRank值（忽略阻尼因子）。

**详细解答**：
1. **构造转移矩阵 $P$**：
   - 网页1有2个外链（到2和3），所以 $P_{21}=0.5$，$P_{31}=0.5$，$P_{11}=0$。
   - 网页2有1个外链（到3），所以 $P_{32}=1$，$P_{12}=0$，$P_{22}=0$。
   - 网页3有1个外链（到1），所以 $P_{13}=1$，$P_{23}=0$，$P_{33}=0$。
   - 转移矩阵 $P = \begin{bmatrix} 0 & 0 & 1 \\ 0.5 & 0 & 0 \\ 0.5 & 1 & 0 \end{bmatrix}$。
2. **初始化**：$r_0 = \begin{bmatrix} 1/3 \\ 1/3 \\ 1/3 \end{bmatrix}$。
3. **第一次迭代**：
   - $r_1 = P r_0 = \begin{bmatrix} 0 & 0 & 1 \\ 0.5 & 0 & 0 \\ 0.5 & 1 & 0 \end{bmatrix} \begin{bmatrix} 1/3 \\ 1/3 \\ 1/3 \end{bmatrix} = \begin{bmatrix} 1/3 \\ 1/6 \\ 1/2 \end{bmatrix}$。
   - 归一化（可选，此处不归一化以观察收敛）。
4. **第二次迭代**：
   - $r_2 = P r_1 = \begin{bmatrix} 0 & 0 & 1 \\ 0.5 & 0 & 0 \\ 0.5 & 1 & 0 \end{bmatrix} \begin{bmatrix} 1/3 \\ 1/6 \\ 1/2 \end{bmatrix} = \begin{bmatrix} 1/2 \\ 1/6 \\ 1/3 \end{bmatrix}$。
5. **继续迭代**：经过多次迭代，$r$ 收敛到 $\begin{bmatrix} 3/7 \\ 2/7 \\ 2/7 \end{bmatrix}$（近似值）。

**答案**：网页1的PageRank值为 $3/7 \approx 0.429$，网页2和3的值均为 $2/7 \approx 0.286$，因此网页1最重要。

---

## 5. 修正例题：幂迭代法

**例题4（修正新增例题3）**：计算矩阵 $A = \begin{bmatrix} 4 & 1 \\ 1 & 3 \end{bmatrix}$ 的最大特征值及其特征向量。

**详细解答**：
1. **初始化**：$x_0 = \begin{bmatrix} 1 \\ 0 \end{bmatrix}$，归一化后仍为 $\begin{bmatrix} 1 \\ 0 \end{bmatrix}$。
2. **第一次迭代**：
   - $x_1 = A x_0 = \begin{bmatrix} 4 & 1 \\ 1 & 3 \end{bmatrix} \begin{bmatrix} 1 \\ 0 \end{bmatrix} = \begin{bmatrix} 4 \\ 1 \end{bmatrix}$。
   - 归一化：$||x_1||_2 = \sqrt{16+1} = \sqrt{17} \approx 4.123$，$x_1 \approx \begin{bmatrix} 0.970 \\ 0.243 \end{bmatrix}$。
3. **第二次迭代**：
   - $x_2 = A x_1 \approx \begin{bmatrix} 4*0.970 + 1*0.243 \\ 1*0.970 + 3*0.243 \end{bmatrix} = \begin{bmatrix} 4.123 \\ 1.699 \end{bmatrix}$。
   - 归一化：$||x_2||_2 \approx 4.456$，$x_2 \approx \begin{bmatrix} 0.925 \\ 0.381 \end{bmatrix}$。
4. **第三次迭代**：
   - $x_3 = A x_2 \approx \begin{bmatrix} 4*0.925 + 1*0.381 \\ 1*0.925 + 3*0.381 \end{bmatrix} = \begin{bmatrix} 4.081 \\ 2.068 \end{bmatrix}$。
   - 归一化：$||x_3||_2 \approx 4.564$，$x_3 \approx \begin{bmatrix} 0.894 \\ 0.453 \end{bmatrix}$。
5. **继续迭代**：向量逐渐收敛到 $\begin{bmatrix} 0.850 \\ 0.527 \end{bmatrix}$，即未归一化形式为 $\begin{bmatrix} 1.615 \\ 1 \end{bmatrix}$ 或近似 $\begin{bmatrix} 1.618 \\ 1 \end{bmatrix}$。
6. **计算特征值**：$\lambda \approx \frac{x_k^T A x_k}{x_k^T x_k}$，取 $x_k \approx \begin{bmatrix} 0.850 \\ 0.527 \end{bmatrix}$：
   - $A x_k \approx \begin{bmatrix} 4*0.850 + 1*0.527 \\ 1*0.850 + 3*0.527 \end{bmatrix} = \begin{bmatrix} 3.927 \\ 2.431 \end{bmatrix}$。
   - $x_k^T A x_k \approx 0.850*3.927 + 0.527*2.431 \approx 4.618$。
   - $x_k^T x_k = 0.850^2 + 0.527^2 \approx 1$。
   - 因此 $\lambda \approx 4.618$。

**答案**：最大特征值约为 $\lambda \approx 4.618$，对应特征向量为 $\begin{bmatrix} 1.618 \\ 1 \end{bmatrix}$（近似值）。

**验证**：特征方程 $\det(A - \lambda I) = (4-\lambda)(3-\lambda) - 1 = \lambda^2 - 7\lambda + 11 = 0$，解得 $\lambda = \frac{7 \pm \sqrt{5}}{2}$，即 $\lambda \approx 4.618$ 或 $\lambda \approx 2.382$，与幂迭代结果一致。

---

## 6. 常见考点与易错点

- **幂迭代法**：初始向量不能与主特征向量正交，否则无法收敛到最大特征值；收敛速度受特征值比值影响。
- **QR算法**：只适合小规模矩阵或对称矩阵，实际应用中需结合位移策略；理解QR分解的正交性和收敛原理。
- **特征值/向量定义与实际意义**：特征值分解在矩阵分析和应用中的重要性。
- **PageRank**：转移矩阵的构造方法，阻尼因子的作用，归一化的必要性。
- **英文术语拼写与公式记忆**：确保考试中能准确写出专业术语和公式。

---

## 7. 英文关键词与表达

- eigenvalue, eigenvector, power iteration, Rayleigh quotient, QR algorithm, diagonalization, convergence, PageRank, principal component analysis (PCA), dominant eigenvalue, damping factor, transition matrix

---

如需更详细例题推导或某一算法的代码实现，可随时补充！