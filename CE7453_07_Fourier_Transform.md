# CE7453 Numerical Algorithms 期末高分超详细攻略（第七章：Fourier Transform）

> **本章内容**：离散傅里叶变换（DFT）、快速傅里叶变换（FFT）、三角插值、JPEG应用  
> **适用对象**：零基础/考前冲刺/快速查漏补缺  
> **内容特色**：详细原理、公式推导、算法流程、例题全解、英文关键词

---

## 1. 基本概念与背景

- **傅里叶变换（Fourier Transform）**：将信号从时域变换到频域，分析信号的频率成分。
- **离散傅里叶变换（DFT）**：对有限长离散信号进行傅里叶变换，常用于数字信号处理、图像压缩等。
- **快速傅里叶变换（FFT）**：高效计算DFT的算法，大幅降低运算量。
- **实际应用**：音频分析、图像处理（如JPEG压缩）、通信、工程仿真等。

---

## 2. 离散傅里叶变换（DFT）

### 2.1 公式与推导

- **DFT公式**：  
  $$
  X_k = \sum_{n=0}^{N-1} x_n e^{-2\pi i kn/N}, \quad k=0,1,...,N-1
  $$
- **IDFT公式**：  
  $$
  x_n = \frac{1}{N} \sum_{k=0}^{N-1} X_k e^{2\pi i kn/N}
  $$
- **性质**：可逆、能量守恒、周期性

### 2.2 计算复杂度

- 直接计算DFT复杂度为 $O(N^2)$，大数据量时效率低

---

## 3. 快速傅里叶变换（FFT）

### 3.1 原理

- **分治思想**：将长度为 $N$ 的DFT分解为两个长度为 $N/2$ 的DFT，递归进行
- **要求**：$N$ 必须为2的幂

### 3.2 算法流程（Cooley-Tukey FFT）

1. 输入序列 $x_0, x_1, ..., x_{N-1}$
2. 若 $N=1$，返回 $x_0$
3. 将序列分为偶数项 $x_{2m}$ 和奇数项 $x_{2m+1}$
4. 分别递归计算偶数和奇数部分的DFT
5. 合并结果：  
   $$
   X_k = E_k + W_N^k O_k \\
   X_{k+N/2} = E_k - W_N^k O_k
   $$
   其中 $E_k$、$O_k$ 分别为偶/奇部分DFT，$W_N^k = e^{-2\pi i k/N}$
6. 复杂度降为 $O(N\log N)$

### 3.3 例题

**例题1**：计算 $x = [1, 0, -1, 0]$ 的DFT

**解答**：
- $N=4$，$X_0 = 1+0+(-1)+0=0$
- $X_1 = 1 - i*(-1) = 1 + i$
- $X_2 = 1+0+(-1)+0=0$
- $X_3 = 1 + i*(-1) = 1 - i$
- 详细可按公式逐项代入

---

## 4. 三角插值与JPEG应用

### 4.1 三角插值（Trigonometric Interpolation）

- 用三角函数（正弦、余弦）拟合周期数据
- DFT本质上就是三角插值系数的计算

### 4.2 JPEG压缩中的DCT

- JPEG采用离散余弦变换（DCT），本质与DFT类似，但只用实数部分
- DCT能高效压缩图像能量，便于去除高频噪声

---

## 5. 常见考点与易错点

- DFT/IDFT公式推导与正负号
- FFT分治递归与合并步骤
- 频域与时域的物理意义
- DCT与DFT的区别
- 英文术语拼写与公式记忆

---

## 6. 英文关键词与表达

- discrete Fourier transform (DFT), fast Fourier transform (FFT), frequency domain, time domain, trigonometric interpolation, discrete cosine transform (DCT), JPEG compression, signal processing, spectrum, complex exponential

---

如需更详细例题推导或某一算法的代码实现，可随时补充！