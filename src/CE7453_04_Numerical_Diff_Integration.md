# CE7453 Numerical Algorithms 期末高分超详细攻略（第四章：Numerical Differentiation & Integration）

> **本章内容**：数值微分、数值积分  
> **适用对象**：零基础/考前冲刺/快速查漏补缺  
> **内容特色**：详细原理、公式推导、算法流程、例题全解、英文关键词

---

## 1. 基本概念与背景

- **数值微分（Numerical Differentiation）**：用离散点近似计算函数的导数，适用于函数表达式复杂或仅有数据点的情况。
- **数值积分（Numerical Integration）**：用有限和近似计算定积分，适用于无法解析积分或仅有数据点的情况。
- **实际应用**：工程仿真、信号处理、物理建模、数据分析等。

---

## 2. 数值微分（Numerical Differentiation）

### 2.1 差分公式（Finite Difference Formulas）

#### 2.1.1 2点前向差分（2-point forward difference）

- **公式**：  
  $$
  f'(x) \approx \frac{f(x+h) - f(x)}{h}
  $$
- **误差**：$O(h)$，截断误差较大

#### 2.1.2 3点中心差分（3-point centered difference）

- **公式**：  
  $$
  f'(x) \approx \frac{f(x+h) - f(x-h)}{2h}
  $$
- **误差**：$O(h^2)$，更精确

#### 2.1.3 5点中心差分（5-point centered difference）

- **公式**：  
  $$
  f'(x) \approx \frac{f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)}{12h}
  $$
- **误差**：$O(h^4)$，高精度

#### 2.1.4 Richardson 外推（Richardson Extrapolation）

- **思想**：用不同步长的近似结果组合，消除主误差项，提高精度。

### 2.2 算法流程（以3点中心差分为例）

1. 选定步长 $h$（过大误差大，过小舍入误差大）
2. 计算 $f(x+h)$ 和 $f(x-h)$
3. 按公式求导数近似值
4. 可用不同 $h$ 验证结果稳定性

### 2.3 典型例题

**例题1**：已知 $f(x)=e^{2x}$，用 $h=0.1$，求 $x=1$ 处的导数近似值（用前向差分法）

**解答**：
- $f(1)=e^{2 \cdot 1}=e^2 \approx 7.3891$
- $f(1+0.1)=e^{2 \cdot 1.1}=e^{2.2} \approx 9.0250$
- $f'(1)\approx \frac{9.0250-7.3891}{0.1}=\frac{1.6359}{0.1} \approx 16.3596$
- 精确值 $f'(1)=2e^2 \approx 14.7782$，误差约为 $1.5814$

**例题2**：已知 $f(x)=e^{2x}$，用 $h=0.1$，求 $x=1$ 处的导数近似值（用中心差分法）

**解答**：
- $f(1+0.1)=e^{2 \cdot 1.1}=e^{2.2} \approx 9.0250$
- $f(1-0.1)=e^{2 \cdot 0.9}=e^{1.8} \approx 6.0496$
- $f'(1)\approx \frac{9.0250-6.0496}{2 \cdot 0.1}=\frac{2.9754}{0.2} \approx 14.8768$
- 精确值 $f'(1)=2e^2 \approx 14.7782$，误差约为 $0.0986$

**例题3**：已知 $f(x)=x^3+2x^2+1$，在 $x=2$ 处用 $h=0.2$，分别用2点前向差分、3点中心差分和5点中心差分计算导数近似值，并与精确值比较。

**解答**：
- 精确导数：$f'(x)=3x^2+4x$，在 $x=2$ 处，$f'(2)=3*(2)^2+4*2=12+8=20$
- **2点前向差分**：
  - $f(2)=2^3+2*2^2+1=8+8+1=17$
  - $f(2+0.2)=f(2.2)=(2.2)^3+2*(2.2)^2+1=10.648+9.68+1=21.328$
  - $f'(2)\approx \frac{21.328-17}{0.2}=\frac{4.328}{0.2}=21.64$
  - 误差：$|21.64-20|=1.64$
- **3点中心差分**：
  - $f(2-0.2)=f(1.8)=(1.8)^3+2*(1.8)^2+1=5.832+6.48+1=13.312$
  - $f(2+0.2)=21.328$ （如上）
  - $f'(2)\approx \frac{21.328-13.312}{2*0.2}=\frac{8.016}{0.4}=20.04$
  - 误差：$|20.04-20|=0.04$
- **5点中心差分**：
  - $f(2-0.4)=f(1.6)=(1.6)^3+2*(1.6)^2+1=4.096+5.12+1=10.216$
  - $f(2-0.2)=13.312$ （如上）
  - $f(2+0.2)=21.328$ （如上）
  - $f(2+0.4)=f(2.4)=(2.4)^3+2*(2.4)^2+1=13.824+11.52+1=26.344$
  - $f'(2)\approx \frac{10.216 - 8*13.312 + 8*21.328 - 26.344}{12*0.2}$
  - $=\frac{10.216 - 106.496 + 170.624 - 26.344}{2.4}=\frac{47.996}{2.4}=19.9983$
  - 误差：$|19.9983-20|=0.0017$
- **比较**：5点中心差分精度最高，误差最小；2点前向差分误差最大。

### 2.4 实际应用与注意事项

- **应用场景**：数值微分常用于速度、加速度计算（如运动学分析），以及优化算法中的梯度计算。
- **注意事项**：
  - 步长 $h$ 选择需平衡截断误差和舍入误差，建议多次尝试不同 $h$ 值。
  - 对于噪声数据，数值微分会放大噪声，需先平滑处理数据。

---

## 3. 数值积分（Numerical Integration）

### 3.1 Newton-Cotes 公式

#### 3.1.1 梯形法（Trapezoidal Rule）

- **公式**：  
  $$
  \int_a^b f(x)dx \approx \frac{h}{2}[f(a) + f(b)]
  $$
  其中 $h=b-a$

- **复合梯形法**：将区间分为 $n$ 段，累加每段梯形面积
  $$
  \int_a^b f(x)dx \approx \frac{h}{2}[f(a) + 2\sum_{k=1}^{n-1}f(a+kh) + f(b)]
  $$

#### 3.1.2 Simpson 法（Simpson's Rule）

- **公式**：  
  $$
  \int_a^b f(x)dx \approx \frac{h}{3}[f(a) + 4f((a+b)/2) + f(b)]
  $$
  $h=(b-a)/2$

- **复合Simpson法**：区间分偶数段，累加每段Simpson近似
  $$
  \int_a^b f(x)dx \approx \frac{h}{3}[f(a) + 4\sum_{k=1,3,...}^{n-1}f(a+kh) + 2\sum_{k=2,4,...}^{n-2}f(a+kh) + f(b)]
  $$

#### 3.1.3 Romberg 积分（Romberg Integration）

- **思想**：用梯形法结果递推外推，提高精度。基于 Richardson 外推，通过多次梯形法计算，逐步消除误差项。

#### 3.1.4 高斯求积（Gaussian Quadrature）

- **思想**：选取最优节点和权重，使多项式积分精度最高。对于 $n$ 点高斯求积，对最高 $2n-1$ 次多项式精确。

### 3.2 算法流程（以复合梯形法为例）

1. 将 $[a,b]$ 分为 $n$ 段，步长 $h=(b-a)/n$
2. 计算端点 $f(a), f(b)$
3. 计算中间点 $f(a+kh)$，$k=1,2,...,n-1$
4. 按公式累加求和

### 3.3 典型例题

**例题1**：用复合梯形法计算 $\int_0^1 e^x dx$，分4段

**解答**：
- $h=0.25$
- $x_0=0, x_1=0.25, x_2=0.5, x_3=0.75, x_4=1$
- $f(x_0)=1, f(x_1)=1.2840, f(x_2)=1.6487, f(x_3)=2.1170, f(x_4)=2.7183$
- $I\approx \frac{0.25}{2}[1+2(1.2840+1.6487+2.1170)+2.7183]=1.7183$
- 精确值 $e-1=1.7183$，误差约为 $0.0000$

**例题2**：用Simpson法计算积分 $\int_0^1 x^2 dx$，分2段

**解答步骤**：
1. $h=0.5$，点：$x_0=0, x_1=0.5, x_2=1$
2. $f(x_0)=0, f(x_1)=0.25, f(x_2)=1$
3. $I \approx \frac{0.5}{3} [0 + 4*0.25 + 1] = \frac{0.5}{3} [0 + 1 + 1] = \frac{0.5}{3} * 2 = \frac{1}{3}$
4. 精确值 $\frac{1}{3}$，误差为0。

**例题3**：用Romberg积分计算 $\int_0^1 \frac{4}{1+x^2} dx$，取初始 $n=1$，进行两次外推。

**解答步骤**：
- 该积分精确值为 $\pi \approx 3.1415926535$，用于验证误差。
- **第一步：梯形法，不同分割**：
  - $n=1$，$h=1$，$T_1(1)=\frac{1}{2}[f(0)+f(1)]=\frac{1}{2}[4+2]=3$
  - $n=2$，$h=0.5$，$T_1(2)=\frac{0.5}{2}[f(0)+2f(0.5)+f(1)]=\frac{0.5}{2}[4+2*3.2+2]=3.2$
  - $n=4$，$h=0.25$，$T_1(4)=\frac{0.25}{2}[f(0)+2(f(0.25)+f(0.5)+f(0.75))+f(1)]$
    - $f(0.25)=\frac{4}{1+0.0625}=3.7647$, $f(0.5)=3.2$, $f(0.75)=\frac{4}{1+0.5625}=2.56$
    - $T_1(4)=\frac{0.25}{2}[4+2*(3.7647+3.2+2.56)+2]=3.1466$
- **第二步：Romberg外推**：
  - 第一次外推：$T_2(2)=\frac{4*T_1(2)-T_1(1)}{4-1}=\frac{4*3.2-3}{3}=3.2667$
  - 第二次外推：$T_2(4)=\frac{4*T_1(4)-T_1(2)}{4-1}=\frac{4*3.1466-3.2}{3}=3.1287$
  - 继续外推：$T_3(4)=\frac{4*T_2(4)-T_2(2)}{4-1}=\frac{4*3.1287-3.2667}{3}=3.1493$
- **结果**：$T_3(4)=3.1493$，与 $\pi$ 较接近，误差约为 $0.0077$。

**例题4**：用三点高斯求积法计算 $\int_{-1}^{1} \sqrt{1-x^2} dx$

**解答步骤**：
1. 三点高斯求积法使用的节点：$x_1 = -\sqrt{3/5}$，$x_2 = 0$，$x_3 = \sqrt{3/5}$
2. 对应的权重：$w_1 = 5/9$，$w_2 = 8/9$，$w_3 = 5/9$
3. 函数值计算：
   - $f(x_1) = f(-\sqrt{3/5}) = \sqrt{1-(\sqrt{3/5})^2} = \sqrt{1-3/5} = \sqrt{2/5} \approx 0.6325$
   - $f(x_2) = f(0) = \sqrt{1-0^2} = 1$
   - $f(x_3) = f(\sqrt{3/5}) = \sqrt{1-(\sqrt{3/5})^2} = \sqrt{1-3/5} = \sqrt{2/5} \approx 0.6325$
4. 求积结果：$I \approx w_1f(x_1) + w_2f(x_2) + w_3f(x_3)$
   - $I \approx (5/9) \cdot 0.6325 + (8/9) \cdot 1 + (5/9) \cdot 0.6325$
   - $I \approx (5/9) \cdot 0.6325 \cdot 2 + (8/9) \cdot 1$
   - $I \approx 0.7028 + 0.8889 \approx 1.5916$
5. 精确值为 $\pi/2 \approx 1.5708$，误差约为 $0.0208$（约1.3%）

### 3.4 实际应用与注意事项

- **应用场景**：数值积分用于计算面积、体积、概率密度积分、物理量累积（如功、能量）。
- **注意事项**：
  - 增加分割数 $n$ 可提高精度，但计算量增加。
  - 对于周期性函数或振荡函数，高斯求积可能更有效。
  - Romberg 积分适合快速收敛，但对函数平滑性要求较高。

---

## 4. 常见考点与易错点

- 步长 $h$ 选取不当导致误差大
- 舍入误差与截断误差权衡
- 复合公式累加时端点权重
- 高斯求积节点和权重记忆
- 英文术语拼写与公式记忆

---

## 5. 英文关键词与表达

- numerical differentiation, finite difference, truncation error, round-off error, Richardson extrapolation, step size
- numerical integration, Newton-Cotes, trapezoidal rule, Simpson's rule, Romberg integration, Gaussian quadrature, composite rule, node, weight

---

## 6. 代码实现示例

下面提供各种数值微分和数值积分方法的Python实现示例。

### 6.1 数值微分代码实现

```python
import numpy as np
import matplotlib.pyplot as plt

def forward_difference(f, x, h=0.01):
    """
    使用前向差分法计算导数
    
    参数:
        f: 函数
        x: 计算导数的点
        h: 步长
    
    返回:
        导数近似值
    """
    return (f(x + h) - f(x)) / h

def central_difference(f, x, h=0.01):
    """
    使用中心差分法计算导数
    
    参数:
        f: 函数
        x: 计算导数的点
        h: 步长
    
    返回:
        导数近似值
    """
    return (f(x + h) - f(x - h)) / (2 * h)

def five_point_difference(f, x, h=0.01):
    """
    使用五点中心差分法计算导数
    
    参数:
        f: 函数
        x: 计算导数的点
        h: 步长
    
    返回:
        导数近似值
    """
    return (f(x - 2*h) - 8*f(x - h) + 8*f(x + h) - f(x + 2*h)) / (12 * h)

def richardson_extrapolation(f, x, h=0.1, k=2):
    """
    使用Richardson外推法提高导数计算精度
    
    参数:
        f: 函数
        x: 计算导数的点
        h: 初始步长
        k: 外推次数
    
    返回:
        改进后的导数近似值
    """
    # 计算中心差分
    D1 = central_difference(f, x, h)
    D2 = central_difference(f, x, h/2)
    
    # 外推公式: D = D2 + (D2 - D1)/(4^k - 1)
    return D2 + (D2 - D1) / (4**k - 1)

# 测试示例
if __name__ == "__main__":
    # 定义测试函数和其精确导数
    f = lambda x: np.exp(2*x)
    df_exact = lambda x: 2 * np.exp(2*x)
    
    # 测试点
    x0 = 1.0
    
    # 计算导数并与精确值比较
    h = 0.1
    fd = forward_difference(f, x0, h)
    cd = central_difference(f, x0, h)
    fpd = five_point_difference(f, x0, h)
    rd = richardson_extrapolation(f, x0, h)
    exact = df_exact(x0)
    
    print(f"函数 f(x) = e^(2x) 在 x = {x0} 处的导数:")
    print(f"精确值: {exact:.6f}")
    print(f"前向差分 (h={h}): {fd:.6f}, 相对误差: {abs(fd-exact)/exact:.6f}")
    print(f"中心差分 (h={h}): {cd:.6f}, 相对误差: {abs(cd-exact)/exact:.6f}")
    print(f"五点差分 (h={h}): {fpd:.6f}, 相对误差: {abs(fpd-exact)/exact:.6f}")
    print(f"Richardson外推: {rd:.6f}, 相对误差: {abs(rd-exact)/exact:.6f}")
    
    # 绘制不同步长下的误差
    h_values = np.logspace(-8, -1, 20)
    fd_errors = [abs(forward_difference(f, x0, h) - exact)/exact for h in h_values]
    cd_errors = [abs(central_difference(f, x0, h) - exact)/exact for h in h_values]
    fpd_errors = [abs(five_point_difference(f, x0, h) - exact)/exact for h in h_values]
    
    plt.figure(figsize=(10, 6))
    plt.loglog(h_values, fd_errors, 'o-', label='前向差分')
    plt.loglog(h_values, cd_errors, 's-', label='中心差分')
    plt.loglog(h_values, fpd_errors, '^-', label='五点差分')
    plt.grid(True)
    plt.xlabel('步长 h')
    plt.ylabel('相对误差')
    plt.title('不同数值微分方法在不同步长下的相对误差')
    plt.legend()
    plt.show()
```

### 6.2 数值积分代码实现

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def trapezoidal_rule(f, a, b, n=100):
    """
    复合梯形法则计算定积分
    
    参数:
        f: 被积函数
        a, b: 积分区间
        n: 分段数
    
    返回:
        积分近似值
    """
    h = (b - a) / n
    x = np.linspace(a, b, n+1)
    y = f(x)
    
    # 复合梯形法则公式
    integral = (h/2) * (y[0] + 2*np.sum(y[1:-1]) + y[-1])
    return integral

def simpson_rule(f, a, b, n=100):
    """
    复合Simpson法则计算定积分
    
    参数:
        f: 被积函数
        a, b: 积分区间
        n: 分段数 (必须为偶数)
    
    返回:
        积分近似值
    """
    if n % 2 != 0:
        n += 1  # 确保n为偶数
    
    h = (b - a) / n
    x = np.linspace(a, b, n+1)
    y = f(x)
    
    # 复合Simpson法则公式
    integral = (h/3) * (y[0] + 4*np.sum(y[1:-1:2]) + 2*np.sum(y[2:-1:2]) + y[-1])
    return integral

def romberg_integration(f, a, b, max_iter=5):
    """
    Romberg积分法
    
    参数:
        f: 被积函数
        a, b: 积分区间
        max_iter: 最大迭代次数
        
    返回:
        积分近似值
    """
    R = np.zeros((max_iter, max_iter))
    
    # 初始计算 - 使用复合梯形法则
    h = b - a
    R[0, 0] = (h/2) * (f(a) + f(b))
    
    for i in range(1, max_iter):
        # 计算下一级梯形法值
        h = h / 2
        n = 2**i
        sum_f = 0
        for k in range(1, n, 2):
            sum_f += f(a + k*h)
        R[i, 0] = R[i-1, 0]/2 + h*sum_f
        
        # Richardson外推
        for j in range(1, i+1):
            R[i, j] = R[i, j-1] + (R[i, j-1] - R[i-1, j-1]) / (4**j - 1)
    
    return R[max_iter-1, max_iter-1]

def gaussian_quadrature(f, a, b, n=5):
    """
    高斯求积法计算定积分
    
    参数:
        f: 被积函数
        a, b: 积分区间
        n: 高斯点数 (最多5点)
        
    返回:
        积分近似值
    """
    # 高斯点和权重 (标准区间 [-1, 1])
    if n == 1:
        x_points = np.array([0])
        weights = np.array([2])
    elif n == 2:
        x_points = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
        weights = np.array([1, 1])
    elif n == 3:
        x_points = np.array([-np.sqrt(3/5), 0, np.sqrt(3/5)])
        weights = np.array([5/9, 8/9, 5/9])
    elif n == 4:
        x0 = np.sqrt((3-2*np.sqrt(6/5))/7)
        x1 = np.sqrt((3+2*np.sqrt(6/5))/7)
        x_points = np.array([-x1, -x0, x0, x1])
        w0 = (18+np.sqrt(30))/36
        w1 = (18-np.sqrt(30))/36
        weights = np.array([w1, w0, w0, w1])
    elif n == 5:
        x0 = 0
        x1 = np.sqrt(5-2*np.sqrt(10/7))/3
        x2 = np.sqrt(5+2*np.sqrt(10/7))/3
        x_points = np.array([-x2, -x1, x0, x1, x2])
        w0 = 128/225
        w1 = (322+13*np.sqrt(70))/900
        w2 = (322-13*np.sqrt(70))/900
        weights = np.array([w2, w1, w0, w1, w2])
    else:
        raise ValueError("目前仅支持1到5个高斯点")
    
    # 转换到积分区间 [a, b]
    scaled_f = lambda x: f((b-a)/2 * x + (a+b)/2) * (b-a)/2
    
    # 计算积分
    integral = np.sum(weights * scaled_f(x_points))
    return integral

# 测试示例
if __name__ == "__main__":
    # 测试函数
    f1 = lambda x: x**2
    f2 = lambda x: np.exp(x)
    f3 = lambda x: 4/(1+x**2)  # 积分结果为π
    f4 = lambda x: np.sqrt(1-x**2)  # 半圆面积，积分结果为π/2
    
    # 积分区间
    a1, b1 = 0, 1  # 用于f1, f2
    a3, b3 = 0, 1  # 用于f3
    a4, b4 = -1, 1  # 用于f4
    
    # 精确结果
    exact1 = 1/3  # ∫(0→1) x^2 dx = 1/3
    exact2 = np.exp(1) - 1  # ∫(0→1) e^x dx = e - 1
    exact3 = np.pi  # ∫(0→1) 4/(1+x^2) dx = π
    exact4 = np.pi/2  # ∫(-1→1) sqrt(1-x^2) dx = π/2
    
    # 不同方法计算结果
    n_intervals = 4  # 分段数
    
    print("函数 f(x) = x^2 在区间 [0, 1] 上的积分:")
    trap1 = trapezoidal_rule(f1, a1, b1, n_intervals)
    simp1 = simpson_rule(f1, a1, b1, n_intervals)
    romb1 = romberg_integration(f1, a1, b1, 3)
    gauss1 = gaussian_quadrature(f1, a1, b1, 3)
    
    print(f"精确值: {exact1}")
    print(f"梯形法则 (n={n_intervals}): {trap1:.6f}, 相对误差: {abs(trap1-exact1)/exact1:.6f}")
    print(f"Simpson法则 (n={n_intervals}): {simp1:.6f}, 相对误差: {abs(simp1-exact1)/exact1:.6f}")
    print(f"Romberg积分: {romb1:.6f}, 相对误差: {abs(romb1-exact1)/exact1:.6f}")
    print(f"高斯求积法 (3点): {gauss1:.6f}, 相对误差: {abs(gauss1-exact1)/exact1:.6f}")
    
    print("\n函数 f(x) = 4/(1+x^2) 在区间 [0, 1] 上的积分 (结果应接近π):")
    trap3 = trapezoidal_rule(f3, a3, b3, n_intervals)
    simp3 = simpson_rule(f3, a3, b3, n_intervals)
    romb3 = romberg_integration(f3, a3, b3, 3)
    gauss3 = gaussian_quadrature(f3, a3, b3, 3)
    
    print(f"精确值: {exact3}")
    print(f"梯形法则 (n={n_intervals}): {trap3:.6f}, 相对误差: {abs(trap3-exact3)/exact3:.6f}")
    print(f"Simpson法则 (n={n_intervals}): {simp3:.6f}, 相对误差: {abs(simp3-exact3)/exact3:.6f}")
    print(f"Romberg积分: {romb3:.6f}, 相对误差: {abs(romb3-exact3)/exact3:.6f}")
    print(f"高斯求积法 (3点): {gauss3:.6f}, 相对误差: {abs(gauss3-exact3)/exact3:.6f}")
    
    # 绘制不同分段数下的误差比较
    n_values = [2, 4, 8, 16, 32, 64, 128]
    
    trap_errors = [abs(trapezoidal_rule(f2, a1, b1, n) - exact2)/exact2 for n in n_values]
    simp_errors = [abs(simpson_rule(f2, a1, b1, n) - exact2)/exact2 for n in n_values]
    
    plt.figure(figsize=(10, 6))
    plt.loglog(n_values, trap_errors, 'o-', label='梯形法则')
    plt.loglog(n_values, simp_errors, 's-', label='Simpson法则')
    plt.grid(True)
    plt.xlabel('分段数 n')
    plt.ylabel('相对误差')
    plt.title('不同数值积分方法在不同分段数下的相对误差 (f(x) = e^x)')
    plt.legend()
    plt.show()
```

### 6.3 应用实例：微分方程数值解法

```python
import numpy as np
import matplotlib.pyplot as plt

def euler_method(f, t0, y0, h, n_steps):
    """
    欧拉法求解常微分方程 dy/dt = f(t, y)
    
    参数:
        f: 函数，表示 dy/dt = f(t, y)
        t0: 初始时间
        y0: 初始值
        h: 步长
        n_steps: 步数
    
    返回:
        t_values: 时间点数组
        y_values: 对应的函数值数组
    """
    t_values = np.zeros(n_steps + 1)
    y_values = np.zeros(n_steps + 1)
    
    t_values[0] = t0
    y_values[0] = y0
    
    for i in range(n_steps):
        t_values[i+1] = t_values[i] + h
        y_values[i+1] = y_values[i] + h * f(t_values[i], y_values[i])
    
    return t_values, y_values

def runge_kutta_4(f, t0, y0, h, n_steps):
    """
    四阶Runge-Kutta方法求解常微分方程 dy/dt = f(t, y)
    
    参数:
        f: 函数，表示 dy/dt = f(t, y)
        t0: 初始时间
        y0: 初始值
        h: 步长
        n_steps: 步数
    
    返回:
        t_values: 时间点数组
        y_values: 对应的函数值数组
    """
    t_values = np.zeros(n_steps + 1)
    y_values = np.zeros(n_steps + 1)
    
    t_values[0] = t0
    y_values[0] = y0
    
    for i in range(n_steps):
        t = t_values[i]
        y = y_values[i]
        
        k1 = f(t, y)
        k2 = f(t + h/2, y + h*k1/2)
        k3 = f(t + h/2, y + h*k2/2)
        k4 = f(t + h, y + h*k3)
        
        t_values[i+1] = t + h
        y_values[i+1] = y + h * (k1 + 2*k2 + 2*k3 + k4) / 6
    
    return t_values, y_values

# 示例：求解微分方程 dy/dt = -y + t (精确解为 y = t - 1 + 2e^(-t))
if __name__ == "__main__":
    # 定义微分方程 dy/dt = f(t, y)
    f = lambda t, y: -y + t
    
    # 定义精确解
    exact_solution = lambda t: t - 1 + 2 * np.exp(-t)
    
    # 初始条件和求解参数
    t0, y0 = 0, 1  # 初始条件
    h = 0.1        # 步长
    n_steps = 50   # 步数
    
    # 使用欧拉法求解
    t_euler, y_euler = euler_method(f, t0, y0, h, n_steps)
    
    # 使用Runge-Kutta法求解
    t_rk4, y_rk4 = runge_kutta_4(f, t0, y0, h, n_steps)
    
    # 计算精确解
    t_exact = np.linspace(t0, t0 + h*n_steps, 200)
    y_exact = exact_solution(t_exact)
    
    # 绘制结果比较
    plt.figure(figsize=(12, 6))
    plt.plot(t_exact, y_exact, 'k-', label='精确解')
    plt.plot(t_euler, y_euler, 'bo-', label='欧拉法')
    plt.plot(t_rk4, y_rk4, 'rs-', label='四阶Runge-Kutta法')
    plt.grid(True)
    plt.xlabel('t')
    plt.ylabel('y')
    plt.title('常微分方程 dy/dt = -y + t 的数值解 (h = ' + str(h) + ')')
    plt.legend()
    
    # 计算并显示误差
    y_exact_at_points = exact_solution(t_euler)
    euler_error = np.abs(y_euler - y_exact_at_points)
    rk4_error = np.abs(y_rk4 - y_exact_at_points)
    
    plt.figure(figsize=(12, 6))
    plt.semilogy(t_euler, euler_error, 'bo-', label='欧拉法误差')
    plt.semilogy(t_euler, rk4_error, 'rs-', label='四阶Runge-Kutta法误差')
    plt.grid(True)
    plt.xlabel('t')
    plt.ylabel('误差 (对数尺度)')
    plt.title('数值解误差比较')
    plt.legend()
    
    plt.show()
```

如需更详细例题推导或某一算法的代码实现，可随时补充！