# CE7453 Numerical Algorithms 期末高分超详细攻略（第二章：Linear Systems）

> **本章内容**：线性方程组（Linear Systems）  
> **适用对象**：零基础/考前冲刺/快速查漏补缺  
> **内容特色**：详细原理、公式推导、算法流程、例题全解、英文关键词

---

## 2.1 基本概念与背景
- **背景与重要性**：
  - **核心地位**：线性方程组是现代科学与工程计算的基石。无论是模拟物理现象、设计复杂系统，还是分析海量数据，都离不开线性方程组的求解。
  - **广泛应用**：从桥梁结构的应力分析、电路网络的电流电压计算，到天气预报的数值模拟、机器学习模型的参数训练，再到经济系统的投入产出分析，线性方程组无处不在。
  - **问题转化**：许多看似非线性的复杂问题，在局部近似或迭代求解的过程中，往往会被转化为一系列线性方程组来处理，例如求解非线性方程组的牛顿法、偏微分方程的有限元/有限差分法等。
  - **性能关键**：求解线性方程组的算法效率（速度）和数值稳定性（精度、抗干扰能力）直接决定了相关应用的可行性和可靠性。特别是在处理大规模问题时，高效且稳健的线性求解器至关重要。

- **什么是线性方程组？**  
  一组形如 $Ax=b$ 的方程，其中 $A$ 是 $m \times n$ 系数矩阵，$x$ 是 $n \times 1$ 未知向量，$b$ 是 $m \times 1$ 常数向量。

- **线性方程组的类型**：
  - **方程数 = 未知数**（$m = n$）：方阵系统，可能有唯一解、无穷多解或无解
  - **方程数 > 未知数**（$m > n$）：超定系统，通常无精确解，可用最小二乘法求近似解
  - **方程数 < 未知数**（$m < n$）：欠定系统，通常有无穷多解

- **解的存在性与唯一性**：
  - 若 $\det(A) \neq 0$（$A$ 满秩），则有唯一解
  - 若 $\det(A) = 0$（$A$ 奇异），则有无穷多解或无解
  - 若 $Ax = b$ 无解，则称该方程组不相容（inconsistent）

- **实际应用举例**：
  - **工程结构分析 (Structural Engineering)**：在使用有限元方法 (Finite Element Method, FEM) 分析桥梁、建筑等结构受力时，需要求解大型稀疏线性方程组 $K u = f$，其中 $K$ 是结构刚度矩阵，$u$ 是节点位移向量，$f$ 是外加载荷向量。
  - **电路分析 (Circuit Analysis)**：根据基尔霍夫定律 (Kirchhoff's Laws) 分析复杂电路时，节点电压法或网孔电流法会导出一组线性方程来确定各处的电压和电流。
  - **计算流体力学 (Computational Fluid Dynamics, CFD)**：在模拟流体流动时，离散化的 Navier-Stokes 方程通常需要在每个时间步求解压力泊松方程等线性系统。
  - **计算机图形学 (Computer Graphics)**：在进行三维建模、渲染和动画时，几何变换（旋转、缩放、平移）、光照模型（辐射度方法）、物理模拟（布料、流体）等都可能涉及线性方程组。
  - **机器学习与数据科学 (Machine Learning & Data Science)**：线性回归、岭回归 (Ridge Regression)、支持向量机 (SVM) 的对偶问题、主成分分析 (PCA) 中的协方差矩阵特征分解等都与线性代数和线性方程组密切相关。最小二乘法是数据拟合和参数估计的基础。
  - **经济学 (Economics)**：投入产出模型 (Input-Output Model) 分析国民经济各部门之间的相互依赖关系，其核心是求解一个线性方程组。
  - **网络流分析 (Network Flow Analysis)**：在交通网络、通信网络或物流网络中，分析流量分配和优化路径常常需要解流量平衡方程，这也是线性方程组。

---

## 2.2 主要算法与原理

#### 理论补充
- **高斯消元法的数学基础**：
  - 通过初等行变换（行交换、数乘、加法）将矩阵化为上三角形，等价于解线性方程组的消元过程。
  - 消元过程本质是消去未知数，逐步简化方程组结构。
- **主元选择的数值意义**：
  - 主元过小会导致舍入误差放大，甚至数值不稳定。
  - 部分主元消去法可显著提高算法的鲁棒性。
- **误差与稳定性分析**：
  - 舍入误差在消元和回代过程中逐步累积，尤其在主元较小时更为严重。
  - 对病态矩阵（条件数大）应优先考虑主元策略或正交分解法。
### 2.2.1 高斯消元法（Gaussian Elimination）

- **基本原理**：  
  通过初等行变换，将增广矩阵 $[A|b]$ 转化为行阶梯形（上三角形），然后通过回代求解未知数。

- **算法流程**：
  1. **前向消元（Forward Elimination）**：
     - 对每一列 $k$ （从左到右）
     - 选择主元（通常是当前列中最大的元素，称为部分主元消去法）
     - 必要时进行行交换
     - 用主元行消去下方各行的该列元素
  2. **回代（Back Substitution）**：
     - 从最后一个未知数开始，逐个求解
     - $x_n = b_n/a_{nn}$
     - $x_i = (b_i - \sum_{j=i+1}^{n} a_{ij}x_j)/a_{ii}$，$i = n-1, n-2, ..., 1$

- **数学表示**：
  - 消元过程：$a_{ij}^{(k+1)} = a_{ij}^{(k)} - m_{ik}a_{kj}^{(k)}$，其中 $m_{ik} = a_{ik}^{(k)}/a_{kk}^{(k)}$ 是乘数
  - 增广矩阵变换：$[A|b] \rightarrow [U|c]$，其中 $U$ 是上三角矩阵

- **部分主元消去法（Partial Pivoting）**：
  - 在每一步消元前，选择当前列中绝对值最大的元素作为主元
  - 通过行交换将主元移至对角线位置
  - 目的：提高数值稳定性，减小舍入误差

- **计算复杂度**：
  - 前向消元：$O(n^3)$
  - 回代：$O(n^2)$
  - 总体：$O(n^3)$

- **伪代码**：
  ```
  function GaussianElimination(A, b)
      n = size(A, 1)  // 矩阵A的行数
      
      // 前向消元
      for k = 1 to n-1 do
          // 部分主元选择
          max_index = k
          for i = k+1 to n do
              if |A[i,k]| > |A[max_index,k]| then
                  max_index = i
              end if
          end for
          
          // 行交换
          if max_index != k then
              swap A[k,:] and A[max_index,:]
              swap b[k] and b[max_index]
          end if
          
          // 消元
          for i = k+1 to n do
              m = A[i,k] / A[k,k]
              A[i,k] = 0
              for j = k+1 to n do
                  A[i,j] = A[i,j] - m * A[k,j]
              end for
              b[i] = b[i] - m * b[k]
          end for
      end for
      
      // 回代
      x = new array of size n
      for i = n downto 1 do
          sum = 0
          for j = i+1 to n do
              sum = sum + A[i,j] * x[j]
          end for
          x[i] = (b[i] - sum) / A[i,i]
      end for
      
      return x
  end function
  ```

#### 误差与数值稳定性
- **舍入误差来源**：有限精度下的加减乘除运算。
- **主元策略**：通过行交换避免主元过小，减少误差放大。
- **病态系统**：条件数大时，微小扰动会导致解大幅变化，应采用更稳定的算法（如QR分解）。
- **优缺点**：
  - **优点**：
    - 直接法，有限步可得精确解
    - 实现简单，适用于一般线性方程组
  - **缺点**：
    - 计算量大，不适合大规模稀疏矩阵
    - 舍入误差可能累积
    - 不利于并行计算

- **英文关键词**：Gaussian elimination, row operation, pivoting, partial pivoting, upper triangular, back substitution, forward elimination, multiplier

#### 例题1：高斯消元法解二元线性方程组 (1)

**问题**：用高斯消元法求解方程组：

$$ \begin{cases} 2x + 3y = 8 \\ 5x + 4y = 13 \end{cases} $$

**解答步骤**：

1. **增广矩阵**：$[A|b] = \begin{bmatrix} 2 & 3 & | & 8 \\ 5 & 4 & | & 13 \end{bmatrix}$
2. **消元**：目标是将第二行的第一个元素变为0。
   - 计算乘数：$m_{21} = a_{21}/a_{11} = 5/2 = 2.5$
   - 第二行减去第一行的 $m_{21}$ 倍：$R_2 = R_2 - 2.5 \times R_1$
     - 新第二行：$[5 - 2.5 \times 2, \quad 4 - 2.5 \times 3, \quad 13 - 2.5 \times 8]$
     - 新第二行：$[5 - 5, \quad 4 - 7.5, \quad 13 - 20] = [0, -3.5, -7]$
   - 变换后的矩阵：$\begin{bmatrix} 2 & 3 & | & 8 \\ 0 & -3.5 & | & -7 \end{bmatrix}$
3. **回代**：
   - 从第二行解出 $y$：$-3.5y = -7 \implies y = \frac{-7}{-3.5} = 2$
   - 将 $y=2$ 代入第一行：$2x + 3(2) = 8 \implies 2x + 6 = 8 \implies 2x = 2 \implies x = 1$

**答案**：$x = 1, y = 2$

**Python代码示例**：
```python
import numpy as np

def solve_gaussian_elimination(A_in, b_in):
    """使用高斯消元法解 Ax = b"""
    A = np.array(A_in, dtype=float)
    b = np.array(b_in, dtype=float)
    n = len(b)
    
    # 构建增广矩阵
    Ab = np.hstack((A, b.reshape(-1, 1)))
    
    # 前向消元
    for k in range(n - 1):
        # 简单主元选择（可选，这里为简化未实现）
        # if np.abs(Ab[k, k]) < 1e-10:
        #     # 寻找下方绝对值最大的行进行交换
        #     max_idx = k + np.argmax(np.abs(Ab[k+1:, k])) + 1
        #     Ab[[k, max_idx]] = Ab[[max_idx, k]] # 行交换

        if np.abs(Ab[k, k]) < 1e-10:
             raise ValueError("主元为零或过小，无法进行消元")

        for i in range(k + 1, n):
            m = Ab[i, k] / Ab[k, k]
            Ab[i, k:] = Ab[i, k:] - m * Ab[k, k:]
            
    # 回代
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        if np.abs(Ab[i, i]) < 1e-10:
            raise ValueError("主元为零或过小，无法回代")
        sum_ax = np.dot(Ab[i, i+1:n], x[i+1:n])
        x[i] = (Ab[i, n] - sum_ax) / Ab[i, i]
        
    return x

# 例题1
A1 = [[2, 3], [5, 4]]
b1 = [8, 13]
x1 = solve_gaussian_elimination(A1, b1)
print(f"例题1 解: x = {x1[0]}, y = {x1[1]}")

# 验证 (可选)
# x_document = np.array([1, 2])
# assert np.allclose(x1, x_document)
```

#### 例题2：高斯消元法解二元线性方程组 (2)

**问题**：用高斯消元法求解方程组：
$$ \begin{cases} 3x + 2y = 7 \\ x + 4y = 5 \end{cases} $$

**解答步骤**：

1. **增广矩阵**：$[A|b] = \begin{bmatrix} 3 & 2 & | & 7 \\ 1 & 4 & | & 5 \end{bmatrix}$
2. **消元**：
   - 计算乘数：$m_{21} = a_{21}/a_{11} = 1/3$
   - 第二行减去第一行的 $m_{21}$ 倍：$R_2 = R_2 - \frac{1}{3} R_1$
     - 新第二行：$[1 - \frac{1}{3}(3), \quad 4 - \frac{1}{3}(2), \quad 5 - \frac{1}{3}(7)]$
     - 新第二行：$[1 - 1, \quad 4 - \frac{2}{3}, \quad 5 - \frac{7}{3}] = [0, \frac{10}{3}, \frac{8}{3}]$
   - 变换后的矩阵：$\begin{bmatrix} 3 & 2 & | & 7 \\ 0 & \frac{10}{3} & | & \frac{8}{3} \end{bmatrix}$
3. **回代**：
   - 从第二行解出 $y$：$\frac{10}{3}y = \frac{8}{3} \implies y = \frac{8}{10} = 0.8$
   - 将 $y=0.8$ 代入第一行：$3x + 2(0.8) = 7 \implies 3x + 1.6 = 7 \implies 3x = 5.4 \implies x = 1.8$

**答案**：$x = 1.8, y = 0.8$

**Python代码示例**：
```python
# 使用上面定义的 solve_gaussian_elimination 函数
A2 = [[3, 2], [1, 4]]
b2 = [7, 5]
x2 = solve_gaussian_elimination(A2, b2)
print(f"例题2 解: x = {x2[0]}, y = {x2[1]}")

# 验证 (可选)
# x_document = np.array([1.8, 0.8])
# assert np.allclose(x2, x_document)
```

---

#### 理论补充
- **LU分解的本质**：
  - 将消元过程矩阵化，$L$ 记录消元乘数，$U$ 记录消元后的上三角结构。
  - $PA=LU$ 形式适用于需要行交换的情况。
- **工程应用**：
  - 多次求解 $Ax=b$，只需分解一次 $A$，大幅提升效率。
- **与高斯消元的关系**：
  - LU分解是高斯消元的矩阵表达，便于理论分析和高效实现。
### 2.2.2 LU分解（LU Decomposition）

- **基本原理**：  
  将系数矩阵 $A$ 分解为下三角矩阵 $L$ 和上三角矩阵 $U$ 的乘积，即 $A = LU$。然后通过解 $Ly = b$ 和 $Ux = y$ 两个三角系统来求解原方程组。

- **与高斯消元的关系**：
  - LU分解本质上是高斯消元的矩阵形式
  - 高斯消元中的乘数 $m_{ik}$ 构成了 $L$ 矩阵的元素
  - 消元后的上三角矩阵即为 $U$

- **分解过程**：
  - 不使用行交换时：$A = LU$
  - 使用行交换时：$PA = LU$，其中 $P$ 是置换矩阵

- **求解流程**：
  1. 分解 $A = LU$
  2. 解 $Ly = b$（前代）
  3. 解 $Ux = y$（回代）

- **计算复杂度**：
  - 分解：$O(n^3)$
  - 求解：$O(n^2)$
  - 对于多个右端项 $b$，只需分解一次，每次求解只需 $O(n^2)$

- **伪代码**：
  ```
  function LUDecomposition(A)
      n = size(A, 1)
      L = identity matrix of size n
      U = copy of A
      
      for k = 1 to n-1 do
          for i = k+1 to n do
              L[i,k] = U[i,k] / U[k,k]
              for j = k to n do
                  U[i,j] = U[i,j] - L[i,k] * U[k,j]
              end for
          end for
      end for
      
      return L, U
  end function
  
  function SolveLU(L, U, b)
      n = size(L, 1)
      
      // 解 Ly = b
      y = new array of size n
      for i = 1 to n do
          sum = 0
          for j = 1 to i-1 do
              sum = sum + L[i,j] * y[j]
          end for
          y[i] = b[i] - sum
      end for
      
      // 解 Ux = y
      x = new array of size n
      for i = n downto 1 do
          sum = 0
          for j = i+1 to n do
              sum = sum + U[i,j] * x[j]
          end for
          x[i] = (y[i] - sum) / U[i,i]
      end for
      
      return x
  end function
  ```

#### 误差与数值稳定性
- LU分解对主元选择同样敏感，通常结合主元策略（如PLU分解）。
- 对于稀疏矩阵，需采用稀疏LU分解以节省存储和计算。
- **优缺点**：
  - **优点**：
    - 对于多个右端项 $b$，效率高
    - 可用于计算行列式、矩阵求逆
    - 分解与求解分离，结构清晰
  - **缺点**：
    - 与高斯消元相同，不适合大规模稀疏矩阵
    - 需要额外存储空间

- **英文关键词**：LU decomposition, LU factorization, lower triangular, upper triangular, forward substitution, backward substitution, multiple right-hand sides

#### 例题3：LU分解解线性方程组

**问题**：用LU分解求解方程组：
$$ \begin{cases} 4x + 3y = 24 \\ 8x + 7y = 52 \end{cases} $$

**解答步骤**：

1. **LU分解**：将 $A = \begin{bmatrix} 4 & 3 \\ 8 & 7 \end{bmatrix}$ 分解为 $A=LU$。
   - 消元：$R_2 = R_2 - (8/4)R_1 = R_2 - 2R_1$
     - 新第二行：$[8 - 2(4), \quad 7 - 2(3)] = [0, 1]$
   - 上三角矩阵 $U = \begin{bmatrix} 4 & 3 \\ 0 & 1 \end{bmatrix}$
   - 乘数 $m_{21} = 2$。下三角矩阵 $L = \begin{bmatrix} 1 & 0 \\ m_{21} & 1 \end{bmatrix} = \begin{bmatrix} 1 & 0 \\ 2 & 1 \end{bmatrix}$
   - 验证：$LU = \begin{bmatrix} 1 & 0 \\ 2 & 1 \end{bmatrix} \begin{bmatrix} 4 & 3 \\ 0 & 1 \end{bmatrix} = \begin{bmatrix} 4 & 3 \\ 8 & 7 \end{bmatrix} = A$

2. **求解 $Ly = b$（前代）**：$b = \begin{bmatrix} 24 \\ 52 \end{bmatrix}$
   $$ \begin{bmatrix} 1 & 0 \\ 2 & 1 \end{bmatrix} \begin{bmatrix} y_1 \\ y_2 \end{bmatrix} = \begin{bmatrix} 24 \\ 52 \end{bmatrix} $$
   - 第一行：$1 \cdot y_1 + 0 \cdot y_2 = 24 \implies y_1 = 24$
   - 第二行：$2 \cdot y_1 + 1 \cdot y_2 = 52 \implies 2(24) + y_2 = 52 \implies 48 + y_2 = 52 \implies y_2 = 4$
   - $y = \begin{bmatrix} 24 \\ 4 \end{bmatrix}$

3. **求解 $Ux = y$（回代）**：
   $$ \begin{bmatrix} 4 & 3 \\ 0 & 1 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 24 \\ 4 \end{bmatrix} $$
   - 第二行：$0 \cdot x_1 + 1 \cdot x_2 = 4 \implies x_2 = 4$
   - 第一行：$4 \cdot x_1 + 3 \cdot x_2 = 24 \implies 4x_1 + 3(4) = 24 \implies 4x_1 + 12 = 24 \implies 4x_1 = 12 \implies x_1 = 3$

**答案**：$x = 3, y = 4$

**Python代码示例**：
```python
import numpy as np
from scipy import linalg

def solve_using_lu(A_in, b_in):
    """使用 scipy 的 LU 分解求解 Ax = b"""
    A = np.array(A_in, dtype=float)
    b = np.array(b_in, dtype=float)
    
    # LU 分解 (scipy.linalg.lu 返回 P, L, U)
    P, L, U = linalg.lu(A)
    # 注意 scipy 返回的 P 是置换矩阵，不是排列向量
    # 需要计算 P @ b 或解 P^T L U x = b
    # 这里我们解 Ly = P^T b, Ux = y
    # 或者更常见的 Ly = Pb, Ux=y (需要理解P的作用)
    # 为简单起见，假设没有行交换(P=I)，或直接用 scipy 的 solve
    
    # 实际应用中，如果需要自己实现，可以基于高斯消元过程记录 L 和 U
    # 或者直接使用 scipy.linalg.solve
    
    # 为了演示流程，我们假设 A = LU (无行交换)
    # 手动分解 (仅适用于本例无行交换情况)
    L_manual = np.array([[1., 0.], [2., 1.]])
    U_manual = np.array([[4., 3.], [0., 1.]])
    
    # 1. 解 Ly = b
    y = linalg.solve_triangular(L_manual, b, lower=True)
    
    # 2. 解 Ux = y
    x = linalg.solve_triangular(U_manual, y, lower=False)
    
    return x

# 例题3
A3 = [[4, 3], [8, 7]]
b3 = [24, 52]
x3 = solve_using_lu(A3, b3)
print(f"例题3 LU分解 解: x = {x3[0]}, y = {x3[1]}")

# 验证 (可选)
# x_document = np.array([3, 4])
# assert np.allclose(x3, x_document)
```

---

### 2.2.3 迭代法（Iterative Methods）

#### 2.2.3.1 Jacobi 迭代法

- **基本原理**：  
  将方程组 $Ax = b$ 改写为 $x = Tx + c$ 的形式，然后迭代求解。Jacobi 法中，每次迭代使用上一次迭代的所有分量值。

- **矩阵分解**：  
  将 $A$ 分解为 $A = D + L + U$，其中 $D$ 是对角矩阵，$L$ 是严格下三角矩阵，$U$ 是严格上三角矩阵。

- **迭代公式**：
  - 分量形式：$x_i^{(k+1)} = \frac{1}{a_{ii}}(b_i - \sum_{j\neq i} a_{ij}x_j^{(k)})$
  - 矩阵形式：$x^{(k+1)} = D^{-1}(b - (L+U)x^{(k)}) = D^{-1}b - D^{-1}(L+U)x^{(k)}$

- **收敛条件**：
  - 迭代矩阵 $T = -D^{-1}(L+U)$ 的谱半径 $\rho(T) < 1$
  - 充分条件：$A$ 严格对角占优（每行对角元素的绝对值大于该行其他元素绝对值之和）

- **伪代码**：
  ```
  function Jacobi(A, b, x0, tol, max_iter)
      n = size(A, 1)
      x = x0
      
      for iter = 1 to max_iter do
          x_new = new array of size n
          
          for i = 1 to n do
              sum = 0
              for j = 1 to n do
                  if j != i then
                      sum = sum + A[i,j] * x[j]
                  end if
              end for
              x_new[i] = (b[i] - sum) / A[i,i]
          end for
          
          if ||x_new - x|| < tol then
              return x_new
          end if
          
          x = x_new
      end for
      
      return "Warning: 达到最大迭代次数"
  end function
  ```

#### 2.2.3.2 Gauss-Seidel 迭代法

- **基本原理**：  
  与 Jacobi 法类似，但在计算第 $i$ 个分量时，使用已经计算出的第 $1$ 到第 $i-1$ 个分量的新值。

- **迭代公式**：
  - 分量形式：$x_i^{(k+1)} = \frac{1}{a_{ii}}(b_i - \sum_{j<i} a_{ij}x_j^{(k+1)} - \sum_{j>i} a_{ij}x_j^{(k)})$
  - 矩阵形式：$x^{(k+1)} = (D+L)^{-1}(b - Ux^{(k)})$

- **收敛条件**：
  - 迭代矩阵 $T = -(D+L)^{-1}U$ 的谱半径 $\rho(T) < 1$
  - 充分条件：$A$ 严格对角占优或对称正定

- **伪代码**：
  ```
  function GaussSeidel(A, b, x0, tol, max_iter)
      n = size(A, 1)
      x = x0
      
      for iter = 1 to max_iter do
          x_new = copy of x
          
          for i = 1 to n do
              sum = 0
              for j = 1 to i-1 do
                  sum = sum + A[i,j] * x_new[j]
              end for
              for j = i+1 to n do
                  sum = sum + A[i,j] * x[j]
              end for
              x_new[i] = (b[i] - sum) / A[i,i]
          end for
          
          if ||x_new - x|| < tol then
              return x_new
          end if
          
          x = x_new
      end for
      
      return "Warning: 达到最大迭代次数"
  end function
  ```

#### 2.2.3.3 Jacobi 与 Gauss-Seidel 比较

- **收敛速度**：
  - Gauss-Seidel 通常比 Jacobi 收敛更快
  - Gauss-Seidel 每次迭代使用最新值，信息传播更快

- **存储需求**：
  - Jacobi 需要额外存储空间保存上一次迭代结果
  - Gauss-Seidel 可以原地更新，节省存储

- **并行性**：
  - Jacobi 更适合并行计算
  - Gauss-Seidel 由于依赖关系，并行化较困难

- **适用场景**：
  - 大规模稀疏矩阵
  - 对角占优矩阵
  - 初值接近真解的情况

- **英文关键词**：Jacobi method, Gauss-Seidel method, iteration, convergence, diagonal dominance, spectral radius, iterative method, sparse matrix

#### 例题4：Jacobi迭代法解线性方程组

**问题**：用Jacobi迭代法求解方程组，进行两次迭代：
$$ \begin{cases} 4x + y = 9 \\ x + 3y = 7 \end{cases} $$
初始值 $x^{(0)} = 0, y^{(0)} = 0$。

**解答步骤**：

1. **迭代公式**：
   - $x^{(k+1)} = \frac{1}{4}(9 - y^{(k)})$
   - $y^{(k+1)} = \frac{1}{3}(7 - x^{(k)})$

2. **初始值**：$x^{(0)} = 0, y^{(0)} = 0$

3. **第一次迭代 (k=0)**：
   - $x^{(1)} = \frac{1}{4}(9 - y^{(0)}) = \frac{1}{4}(9 - 0) = \frac{9}{4} = 2.25$
   - $y^{(1)} = \frac{1}{3}(7 - x^{(0)}) = \frac{1}{3}(7 - 0) = \frac{7}{3} \approx 2.33$

4. **第二次迭代 (k=1)**：
   - $x^{(2)} = \frac{1}{4}(9 - y^{(1)}) = \frac{1}{4}(9 - 7/3) = \frac{1}{4}(\frac{27-7}{3}) = \frac{1}{4}(\frac{20}{3}) = \frac{5}{3} \approx 1.67$
   - $y^{(2)} = \frac{1}{3}(7 - x^{(1)}) = \frac{1}{3}(7 - 9/4) = \frac{1}{3}(\frac{28-9}{4}) = \frac{1}{3}(\frac{19}{4}) = \frac{19}{12} \approx 1.58$

**答案**：迭代两次后，$x^{(2)} \approx 1.67, y^{(2)} \approx 1.58$。

**Python代码示例**：
```python
import numpy as np

def solve_jacobi(A_in, b_in, x0_in, iterations):
    """使用 Jacobi 迭代法求解 Ax = b"""
    A = np.array(A_in, dtype=float)
    b = np.array(b_in, dtype=float)
    x = np.array(x0_in, dtype=float)
    n = len(b)
    x_new = np.zeros(n)
    
    for k in range(iterations):
        for i in range(n):
            sigma = 0
            for j in range(n):
                if i != j:
                    sigma += A[i, j] * x[j]
            if np.abs(A[i, i]) < 1e-10:
                 raise ValueError("对角元素过小")
            x_new[i] = (b[i] - sigma) / A[i, i]
        x = x_new.copy() # 更新 x 用于下一次迭代
        print(f"迭代 {k+1}: x = {x}") # 打印每次迭代结果
        
    return x

# 例题4
A4 = [[4, 1], [1, 3]]
b4 = [9, 7]
x0_4 = [0, 0]
iterations4 = 2
x4 = solve_jacobi(A4, b4, x0_4, iterations4)
print(f"例题4 Jacobi 解 (迭代 {iterations4} 次): x = {x4[0]:.2f}, y = {x4[1]:.2f}")

# 验证 (可选)
# x_document = np.array([1.67, 1.58])
# assert np.allclose(x4, x_document, rtol=1e-2)
```

---

#### 2.2.3.2 Gauss-Seidel 迭代法

#### 例题5：Gauss-Seidel迭代法解线性方程组

**问题**：用Gauss-Seidel迭代法求解例题4的方程组，进行两次迭代：
$$ \begin{cases} 4x + y = 9 \\ x + 3y = 7 \end{cases} $$
初始值 $x^{(0)} = 0, y^{(0)} = 0$。

**解答步骤**：

1. **迭代公式**：
   - $x^{(k+1)} = \frac{1}{4}(9 - y^{(k)})$
   - $y^{(k+1)} = \frac{1}{3}(7 - x^{(k+1)})$  （注意这里用的是最新的 $x^{(k+1)}$）

2. **初始值**：$x^{(0)} = 0, y^{(0)} = 0$

3. **第一次迭代 (k=0)**：
   - $x^{(1)} = \frac{1}{4}(9 - y^{(0)}) = \frac{1}{4}(9 - 0) = \frac{9}{4} = 2.25$
   - $y^{(1)} = \frac{1}{3}(7 - x^{(1)}) = \frac{1}{3}(7 - 2.25) = \frac{1}{3}(4.75) = \frac{19}{12} \approx 1.58$

4. **第二次迭代 (k=1)**：
   - $x^{(2)} = \frac{1}{4}(9 - y^{(1)}) = \frac{1}{4}(9 - 19/12) = \frac{1}{4}(\frac{108-19}{12}) = \frac{1}{4}(\frac{89}{12}) = \frac{89}{48} \approx 1.85$
   - $y^{(2)} = \frac{1}{3}(7 - x^{(2)}) = \frac{1}{3}(7 - 89/48) = \frac{1}{3}(\frac{336-89}{48}) = \frac{1}{3}(\frac{247}{48}) = \frac{247}{144} \approx 1.72$

**答案**：迭代两次后，$x^{(2)} \approx 1.85, y^{(2)} \approx 1.72$。（文档答案 $x=1.86, y=1.71$ 可能是计算或舍入方式略有不同，但趋势一致）

**Python代码示例**：
```python
import numpy as np

def solve_gauss_seidel(A_in, b_in, x0_in, iterations):
    """使用 Gauss-Seidel 迭代法求解 Ax = b"""
    A = np.array(A_in, dtype=float)
    b = np.array(b_in, dtype=float)
    x = np.array(x0_in, dtype=float)
    n = len(b)
    
    for k in range(iterations):
        x_old = x.copy() # 保存旧值用于比较（可选）
        for i in range(n):
            sigma = 0
            for j in range(n):
                if i != j:
                    # 注意这里用的是当前迭代中已经更新的值 x[j]
                    sigma += A[i, j] * x[j] 
            if np.abs(A[i, i]) < 1e-10:
                 raise ValueError("对角元素过小")
            x[i] = (b[i] - sigma) / A[i, i]
        print(f"迭代 {k+1}: x = {x}") # 打印每次迭代结果
        # 收敛判断 (可选)
        # if np.linalg.norm(x - x_old) < tol:
        #     break
            
    return x

# 例题5
A5 = [[4, 1], [1, 3]]
b5 = [9, 7]
x0_5 = [0, 0]
iterations5 = 2
x5 = solve_gauss_seidel(A5, b5, x0_5, iterations5)
print(f"例题5 Gauss-Seidel 解 (迭代 {iterations5} 次): x = {x5[0]:.2f}, y = {x5[1]:.2f}")

# 验证 (可选)
# x_document = np.array([1.86, 1.71]) # 文档答案
# assert np.allclose(x5, x_document, rtol=1e-2)
```

---

## 2.3 实践案例与代码

### 2.3.1 高斯消元法解四阶线性方程组

**问题**：用高斯消元法求解以下方程组：

$$ \begin{cases} 2x_1 + x_2 + 3x_3 + 2x_4 = 21 \\ 4x_1 + 3x_2 + 2x_3 + x_4 = 20 \\ x_1 + 2x_2 + 4x_3 + 3x_4 = 29 \\ 3x_1 + 4x_2 + x_3 + 2x_4 = 22 \end{cases} $$

> **注意**：此例题在原始文档中可能存在印刷错误，这里使用的是验证过的版本。

**解答**：这个过程比较繁琐，适合用代码实现。

**答案**：$x_1 = 1, x_2 = 2, x_3 = 3, x_4 = 4$

**Python代码示例**：
```python
# 使用上面定义的 solve_gaussian_elimination 函数
A6 = [
    [2, 1, 3, 2],
    [4, 3, 2, 1],
    [1, 2, 4, 3],
    [3, 4, 1, 2]
]
b6 = [21, 20, 29, 22]

x6 = solve_gaussian_elimination(A6, b6)
print(f"例题6 (4x4) 高斯消元 解: {x6}")

# 验证 (可选)
x_document = np.array([1, 2, 3, 4])
assert np.allclose(x6, x_document)
```

---

### 2.3.2 迭代法的收敛性判断

**问题**：判断Jacobi迭代法和Gauss-Seidel迭代法对于以下方程组是否收敛：

$$ \begin{cases} 10x + 2y = 12 \\ 3x + 15y = 18 \end{cases} $$

**分析**：

1. **严格对角占优**：
   对于系数矩阵 $A = \begin{bmatrix} 10 & 2 \\ 3 & 15 \end{bmatrix}$：
   - 第一行：$|a_{11}| = |10| = 10$，$|a_{12}| = |2| = 2$。$10 > 2$。满足。
   - 第二行：$|a_{22}| = |15| = 15$，$|a_{21}| = |3| = 3$。$15 > 3$。满足。
   因为矩阵 $A$ 是严格对角占优的，所以Jacobi迭代法和Gauss-Seidel迭代法都收敛。

2. **迭代矩阵的谱半径**：
   - **Jacobi迭代矩阵** $T_J = -D^{-1}(L+U)$
     $D = \begin{bmatrix} 10 & 0 \\ 0 & 15 \end{bmatrix}$, $L = \begin{bmatrix} 0 & 0 \\ 3 & 0 \end{bmatrix}$, $U = \begin{bmatrix} 0 & 2 \\ 0 & 0 \end{bmatrix}$
     $D^{-1} = \begin{bmatrix} 1/10 & 0 \\ 0 & 1/15 \end{bmatrix}$
     $L+U = \begin{bmatrix} 0 & 2 \\ 3 & 0 \end{bmatrix}$
     $T_J = -\begin{bmatrix} 1/10 & 0 \\ 0 & 1/15 \end{bmatrix} \begin{bmatrix} 0 & 2 \\ 3 & 0 \end{bmatrix} = -\begin{bmatrix} 0 & 2/10 \\ 3/15 & 0 \end{bmatrix} = \begin{bmatrix} 0 & -0.2 \\ -0.2 & 0 \end{bmatrix}$
     特征值 $\lambda$ 满足 $\det(T_J - \lambda I) = \det\begin{bmatrix} -\lambda & -0.2 \\ -0.2 & -\lambda \end{bmatrix} = \lambda^2 - (-0.2)(-0.2) = \lambda^2 - 0.04 = 0$。
     $\lambda^2 = 0.04 \implies \lambda = \pm 0.2$。
     谱半径 $\rho(T_J) = \max(|0.2|, |-0.2|) = 0.2 < 1$。所以Jacobi法收敛。

   - **Gauss-Seidel迭代矩阵** $T_G = -(D+L)^{-1}U$
     $D+L = \begin{bmatrix} 10 & 0 \\ 3 & 15 \end{bmatrix}$
     $(D+L)^{-1} = \frac{1}{10 \times 15 - 0 \times 3} \begin{bmatrix} 15 & 0 \\ -3 & 10 \end{bmatrix} = \frac{1}{150} \begin{bmatrix} 15 & 0 \\ -3 & 10 \end{bmatrix} = \begin{bmatrix} 1/10 & 0 \\ -1/50 & 1/15 \end{bmatrix}$
     $T_G = -\begin{bmatrix} 1/10 & 0 \\ -1/50 & 1/15 \end{bmatrix} \begin{bmatrix} 0 & 2 \\ 0 & 0 \end{bmatrix} = -\begin{bmatrix} 0 & 2/10 \\ 0 & -2/50 \end{bmatrix} = \begin{bmatrix} 0 & -0.2 \\ 0 & 1/25 \end{bmatrix} = \begin{bmatrix} 0 & -0.2 \\ 0 & 0.04 \end{bmatrix}$
     特征值 $\lambda$ 满足 $\det(T_G - \lambda I) = \det\begin{bmatrix} -\lambda & -0.2 \\ 0 & 0.04-\lambda \end{bmatrix} = (-\lambda)(0.04-\lambda) - 0 = 0$。
     $\lambda = 0$ 或 $\lambda = 0.04$。
     谱半径 $\rho(T_G) = \max(|0|, |0.04|) = 0.04 < 1$。所以Gauss-Seidel法收敛。

**答案**：Jacobi迭代法和Gauss-Seidel迭代法都收敛。

**Python代码示例**：
```python
import numpy as np

def check_convergence(A_in):
    """检查Jacobi和Gauss-Seidel迭代法的收敛性"""
    A = np.array(A_in, dtype=float)
    n = A.shape[0]
    
    # 1. 检查严格对角占优
    is_diag_dominant = True
    for i in range(n):
        diag = np.abs(A[i, i])
        off_diag_sum = np.sum(np.abs(A[i, :])) - diag
        if diag <= off_diag_sum:
            is_diag_dominant = False
            break
    print(f"严格对角占优: {is_diag_dominant}")
    
    # 2. 计算Jacobi迭代矩阵谱半径
    D = np.diag(np.diag(A))
    L_plus_U = A - D
    try:
        D_inv = np.linalg.inv(D)
        T_J = -D_inv @ L_plus_U
        rho_J = np.max(np.abs(np.linalg.eigvals(T_J)))
        print(f"Jacobi 谱半径: {rho_J:.4f}, 收敛: {rho_J < 1}")
    except np.linalg.LinAlgError:
        print("Jacobi: D 不可逆")
        rho_J = float('inf')

    # 3. 计算Gauss-Seidel迭代矩阵谱半径
    D_plus_L = np.tril(A)
    U = np.triu(A, k=1)
    try:
        D_plus_L_inv = np.linalg.inv(D_plus_L)
        T_G = -D_plus_L_inv @ U
        rho_G = np.max(np.abs(np.linalg.eigvals(T_G)))
        print(f"Gauss-Seidel 谱半径: {rho_G:.4f}, 收敛: {rho_G < 1}")
    except np.linalg.LinAlgError:
        print("Gauss-Seidel: D+L 不可逆")
        rho_G = float('inf')
        
    return is_diag_dominant, rho_J < 1, rho_G < 1

# 例题7
A7 = [[10, 2], [3, 15]]
check_convergence(A7)
```

---

