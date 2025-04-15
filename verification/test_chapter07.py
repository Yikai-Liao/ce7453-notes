import numpy as np
import pytest

def is_close(a, b, rtol=1e-2, atol=1e-8):
    """检查两个数值是否足够接近，考虑浮点数精度问题"""
    if isinstance(a, complex) or isinstance(b, complex):
        return np.abs(a - b) <= (atol + rtol * np.abs(b))
    return np.abs(a - b) <= (atol + rtol * np.abs(b))

def dft(x):
    """
    计算离散傅里叶变换（DFT）
    参数：
        x: 输入信号（列表或数组）
    返回：
        X: DFT结果（复数数组）
    """
    N = len(x)
    X = np.zeros(N, dtype=complex)
    
    for k in range(N):
        for n in range(N):
            X[k] += x[n] * np.exp(-2j * np.pi * k * n / N)
    
    return X

def test_dft_example1():
    """
    验证例题1：计算 x = [1, 0, -1, 0] 的 DFT
    文档中的答案：X = [0, 2, 0, 2]
    """
    # 输入信号
    x = np.array([1, 0, -1, 0])
    
    # 使用自定义DFT函数计算
    X_custom = dft(x)
    
    # 使用NumPy库计算DFT进行验证
    X_numpy = np.fft.fft(x)
    
    # 期望结果
    expected_X = np.array([0, 2, 0, 2], dtype=complex)
    
    # 打印结果进行比较
    print("自定义DFT计算结果:", X_custom)
    print("NumPy FFT计算结果:", X_numpy)
    print("期望结果:", expected_X)
    
    # 验证结果是否与期望值一致
    for i in range(len(expected_X)):
        assert is_close(X_custom[i], expected_X[i])
        assert is_close(X_numpy[i], expected_X[i])

def test_dft_example3():
    """
    验证例题3：计算信号 x = [1, 1, 1, 1] 的 DFT
    文档中的答案：X = [4, 0, 0, 0]
    """
    # 输入信号
    x = np.array([1, 1, 1, 1])
    
    # 使用自定义DFT函数计算
    X_custom = dft(x)
    
    # 使用NumPy库计算DFT进行验证
    X_numpy = np.fft.fft(x)
    
    # 期望结果
    expected_X = np.array([4, 0, 0, 0], dtype=complex)
    
    # 打印结果进行比较
    print("自定义DFT计算结果:", X_custom)
    print("NumPy FFT计算结果:", X_numpy)
    print("期望结果:", expected_X)
    
    # 验证结果是否与期望值一致
    for i in range(len(expected_X)):
        assert is_close(X_custom[i], expected_X[i])
        assert is_close(X_numpy[i], expected_X[i])

def test_fft_example2():
    """
    验证例题2：使用 FFT 计算 x = [1, 2, 3, 4, 5, 6, 7, 8] 的 DFT
    文档中的答案：X = [36, -4+9.656i, -4+4i, -4+1.656i, -4, -4-1.656i, -4-4i, -4-9.656i]
    """
    # 输入信号
    x = np.array([1, 2, 3, 4, 5, 6, 7, 8])
    
    # 使用NumPy FFT库计算
    X_numpy = np.fft.fft(x)
    
    # 期望结果（根据文档中提供的答案）
    expected_X = np.array([
        36, 
        -4 + 9.656j, 
        -4 + 4j, 
        -4 + 1.656j, 
        -4, 
        -4 - 1.656j, 
        -4 - 4j, 
        -4 - 9.656j
    ])
    
    # 打印结果进行比较
    print("NumPy FFT计算结果:", X_numpy)
    print("期望结果:", expected_X)
    
    # 验证结果是否与期望值一致（使用较宽松的误差范围）
    for i in range(len(expected_X)):
        assert is_close(X_numpy[i], expected_X[i], rtol=1e-1)

def idft(X):
    """
    计算逆离散傅里叶变换（IDFT）
    参数：
        X: 频域信号（列表或数组）
    返回：
        x: IDFT结果（复数数组）
    """
    N = len(X)
    x = np.zeros(N, dtype=complex)
    
    for n in range(N):
        for k in range(N):
            x[n] += X[k] * np.exp(2j * np.pi * k * n / N)
        x[n] /= N
    
    return x

def test_dft_idft_property():
    """
    验证DFT和IDFT的可逆性质
    """
    # 测试信号
    signals = [
        np.array([1, 0, -1, 0]),  # 例题1中的信号
        np.array([1, 1, 1, 1]),   # 例题3中的信号
        np.array([1, 2, 3, 4, 5, 6, 7, 8])  # 例题2中的信号
    ]
    
    for x in signals:
        # 应用DFT
        X = dft(x)
        # 应用IDFT
        x_recovered = idft(X)
        
        # 验证恢复的信号是否与原信号一致
        for i in range(len(x)):
            assert is_close(x[i], x_recovered[i])
        
        # 使用NumPy库进行验证
        X_numpy = np.fft.fft(x)
        x_recovered_numpy = np.fft.ifft(X_numpy)
        
        for i in range(len(x)):
            assert is_close(x[i], x_recovered_numpy[i])

if __name__ == "__main__":
    # 直接调用测试函数
    print("测试DFT例题1...")
    test_dft_example1()
    print("\n测试DFT例题3...")
    test_dft_example3()
    print("\n测试FFT例题2...")
    test_fft_example2()
    print("\n测试DFT和IDFT的可逆性质...")
    test_dft_idft_property()
    print("\n所有测试通过!") 