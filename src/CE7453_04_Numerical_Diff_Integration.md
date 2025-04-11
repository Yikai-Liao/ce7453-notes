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

**例题1**：已知 $f(x)=\sin(x)$，用 $h=0.1$，求 $x=1$ 处的导数近似值（用3点中心差分）

**解答**：
- $f(1+0.1)=\sin(1.1)\approx 0.8912$
- $f(1-0.1)=\sin(0.9)\approx 0.7833$
- $f'(1)\approx \frac{0.8912-0.7833}{0.2}=0.5395$
- 精确值 $f'(1)=\cos(1)\approx 0.5403$，误差很小

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

#### 3.1.2 Simpson 法（Simpson's Rule）

- **公式**：  
  $$
  \int_a^b f(x)dx \approx \frac{h}{3}[f(a) + 4f((a+b)/2) + f(b)]
  $$
  $h=(b-a)/2$

- **复合Simpson法**：区间分偶数段，累加每段Simpson近似

#### 3.1.3 Romberg 积分（Romberg Integration）

- **思想**：用梯形法结果递推外推，提高精度

#### 3.1.4 高斯求积（Gaussian Quadrature）

- **思想**：选取最优节点和权重，使多项式积分精度最高

### 3.2 算法流程（以复合梯形法为例）

1. 将 $[a,b]$ 分为 $n$ 段，步长 $h=(b-a)/n$
2. 计算端点 $f(a), f(b)$
3. 计算中间点 $f(a+kh)$，$k=1,2,...,n-1$
4. 按公式累加求和

### 3.3 典型例题

**例题2**：用复合梯形法计算 $\int_0^1 e^x dx$，分4段

**解答**：
- $h=0.25$
- $x_0=0, x_1=0.25, x_2=0.5, x_3=0.75, x_4=1$
- $f(x_0)=1, f(x_1)=1.2840, f(x_2)=1.6487, f(x_3)=2.1170, f(x_4)=2.7183$
- $I\approx \frac{0.25}{2}[1+2(1.2840+1.6487+2.1170)+2.7183]=1.7183$
- 精确值 $e-1\approx 1.7183$

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

如需更详细例题推导或某一算法的代码实现，可随时补充！