import numpy as np
import pytest
from scipy import integrate

def is_close(a, b, rtol=1e-5, atol=1e-8):
    """
    检查两个数值是否足够接近，考虑浮点数精度问题
    """
    return np.abs(a - b) <= (atol + rtol * np.abs(b))

def test_forward_difference():
    """
    验证前向差分法计算导数的例题
    f(x) = e^(2x), x = 1, h = 0.1
    文档中的答案：f'(1) ≈ 16.3596
    """
    # 定义函数
    f = lambda x: np.exp(2*x)
    x0 = 1
    h = 0.1
    
    # 前向差分计算导数
    df_forward = (f(x0 + h) - f(x0)) / h
    
    # 文档中给出的答案
    expected_derivative = 16.3596
    
    # 验证结果是否与文档答案一致
    assert is_close(df_forward, expected_derivative, rtol=1e-4)
    
    # 计算精确导数进行比较
    exact_derivative = 2 * np.exp(2*x0)  # f'(x) = 2e^(2x)
    
    # 验证文档答案相对于精确值的误差
    error = abs(expected_derivative - exact_derivative) / exact_derivative
    assert error < 0.12  # 相对误差应小于12%

def test_central_difference():
    """
    验证中心差分法计算导数的例题
    f(x) = e^(2x), x = 1, h = 0.1
    文档中的答案：f'(1) ≈ 14.8768
    """
    # 定义函数
    f = lambda x: np.exp(2*x)
    x0 = 1
    h = 0.1
    
    # 中心差分计算导数
    df_central = (f(x0 + h) - f(x0 - h)) / (2*h)
    
    # 文档中给出的答案
    expected_derivative = 14.8768
    
    # 验证结果是否与文档答案一致
    assert is_close(df_central, expected_derivative, rtol=1e-4)
    
    # 计算精确导数进行比较
    exact_derivative = 2 * np.exp(2*x0)  # f'(x) = 2e^(2x)
    
    # 验证中心差分相对于前向差分的精度提高
    forward_diff = (f(x0 + h) - f(x0)) / h
    central_error = abs(df_central - exact_derivative) / exact_derivative
    forward_error = abs(forward_diff - exact_derivative) / exact_derivative
    assert central_error < forward_error

def test_richardson_extrapolation():
    """
    验证Richardson外推法计算导数的例题
    f(x) = sin(x), x = π/4, h1 = 0.1, h2 = 0.05
    文档中的答案：f'(π/4) ≈ 0.7071
    """
    # 定义函数
    f = lambda x: np.sin(x)
    x0 = np.pi/4
    h1 = 0.1
    h2 = 0.05
    
    # 使用两个步长计算中心差分
    df_h1 = (f(x0 + h1) - f(x0 - h1)) / (2*h1)
    df_h2 = (f(x0 + h2) - f(x0 - h2)) / (2*h2)
    
    # Richardson外推
    df_richardson = df_h2 + (df_h2 - df_h1) / (2**2 - 1)
    
    # 文档中给出的答案
    expected_derivative = 0.7071
    
    # 验证结果是否与文档答案一致
    assert is_close(df_richardson, expected_derivative, rtol=1e-4)
    
    # 计算精确导数进行比较 (sin(x)的导数是cos(x))
    exact_derivative = np.cos(x0)
    
    # 验证外推法相对于简单中心差分的精度提高
    assert abs(df_richardson - exact_derivative) < abs(df_h2 - exact_derivative)

def test_trapezoidal_rule():
    """
    验证梯形法则计算定积分的例题
    ∫(0→1) x^2 dx, n = 4
    文档中的答案：∫(0→1) x^2 dx ≈ 0.34375
    """
    # 定义函数
    f = lambda x: x**2
    a, b = 0, 1
    n = 4
    
    # 计算步长
    h = (b - a) / n
    
    # 使用梯形法则
    x = np.linspace(a, b, n+1)
    y = f(x)
    integral_trap = (h/2) * (y[0] + 2*np.sum(y[1:-1]) + y[-1])
    
    # 文档中给出的答案
    expected_integral = 0.34375
    
    # 验证结果是否与文档答案一致
    assert is_close(integral_trap, expected_integral, rtol=1e-4)
    
    # 计算精确积分进行比较
    exact_integral = 1/3  # ∫(0→1) x^2 dx = 1/3
    
    # 验证误差是否在合理范围内
    error = abs(integral_trap - exact_integral)
    assert error < 0.011  # 允许更大的误差范围

def test_simpson_rule():
    """
    验证Simpson法则计算定积分的例题
    ∫(0→1) x^2 dx, n = 4
    文档中的答案：∫(0→1) x^2 dx ≈ 0.3333
    """
    # 定义函数
    f = lambda x: x**2
    a, b = 0, 1
    n = 4  # 必须是偶数
    
    # 计算步长
    h = (b - a) / n
    
    # 使用Simpson法则
    x = np.linspace(a, b, n+1)
    y = f(x)
    integral_simpson = (h/3) * (y[0] + 4*np.sum(y[1:-1:2]) + 2*np.sum(y[2:-1:2]) + y[-1])
    
    # 文档中给出的答案
    expected_integral = 0.3333
    
    # 验证结果是否与文档答案一致
    assert is_close(integral_simpson, expected_integral, rtol=1e-4)
    
    # 计算精确积分进行比较
    exact_integral = 1/3  # ∫(0→1) x^2 dx = 1/3
    
    # 验证Simpson法则相对于梯形法则的精度提高
    trap_integral = (h/2) * (y[0] + 2*np.sum(y[1:-1]) + y[-1])
    simpson_error = abs(integral_simpson - exact_integral)
    trap_error = abs(trap_integral - exact_integral)
    assert simpson_error < trap_error

def test_gaussian_quadrature():
    """
    验证高斯求积法计算定积分的例题
    ∫(-1→1) sqrt(1-x^2) dx, n=3
    文档中的答案：∫(-1→1) sqrt(1-x^2) dx ≈ 1.5916 (π/2 ≈ 1.5708)
    """
    # 定义函数 (f(x) = sqrt(1-x^2))
    f = lambda x: np.sqrt(1-x**2)
    
    # 三点高斯求积法
    # 高斯点和权重 (3点公式)
    x_points = np.array([-np.sqrt(3/5), 0, np.sqrt(3/5)])
    weights = np.array([5/9, 8/9, 5/9])
    
    # 计算积分
    integral_gauss = np.sum(weights * f(x_points))
    
    # 文档中给出的答案
    expected_integral = 1.5916
    
    # 验证结果是否与文档答案一致
    assert is_close(integral_gauss, expected_integral, rtol=1e-3)
    
    # 计算精确积分进行比较 (这是半圆面积，值为π/2)
    exact_integral = np.pi/2
    
    # 验证误差是否在合理范围内
    error = abs(integral_gauss - exact_integral) / exact_integral
    assert error < 0.02

if __name__ == "__main__":
    # 运行所有测试
    pytest.main(["-v", __file__]) 