import numpy as np
import pytest
from scipy import optimize

# 定义辅助函数
def is_close(a, b, rtol=1e-5, atol=1e-8):
    """
    检查两个数值是否足够接近，考虑浮点数精度问题
    """
    return np.abs(a - b) <= (atol + rtol * np.abs(b))

def is_array_close(a, b, rtol=1e-5, atol=1e-8):
    """
    检查两个数组是否足够接近，考虑浮点数精度问题
    """
    if np.array(a).shape != np.array(b).shape:
        return False
    return np.allclose(a, b, rtol=rtol, atol=atol)

def relative_error(computed, exact):
    """
    计算相对误差
    """
    if np.isscalar(exact) and exact == 0:
        return np.abs(computed)
    elif hasattr(exact, '__len__') and np.allclose(exact, 0):
        return np.linalg.norm(computed)
    else:
        if np.isscalar(exact):
            return np.abs(computed - exact) / np.abs(exact)
        else:
            return np.linalg.norm(computed - exact) / np.linalg.norm(exact)

def test_example1_bisection():
    """
    验证例题1：用二分法求 f(x)=x^3+x-1 在 [0,1] 内的根，精度 10^-2
    文档中的答案：x ≈ 0.6796875
    """
    f = lambda x: x**3 + x - 1
    
    # 二分法实现
    def bisection_method(f, a, b, tol=1e-10, max_iter=100):
        if f(a) * f(b) >= 0:
            raise ValueError("区间端点函数值必须异号")
        
        iter_count = 0
        while (b - a) / 2 > tol and iter_count < max_iter:
            c = (a + b) / 2
            if f(c) == 0:
                return c
            elif f(a) * f(c) < 0:
                b = c
            else:
                a = c
            iter_count += 1
        
        return (a + b) / 2
    
    # 用自己实现的二分法计算
    root_bisection = bisection_method(f, 0, 1, tol=1e-2)
    # 用scipy的方法计算精确值
    root_exact = optimize.root_scalar(f, bracket=[0, 1], method='brentq').root
    
    # 文档中给出的答案
    root_document = 0.6796875
    
    # 验证自己实现的方法是否与文档答案一致
    assert is_close(root_bisection, root_document, rtol=1e-5)
    # 验证文档答案是否足够接近精确值
    assert is_close(root_document, root_exact, rtol=2e-2)
    # 验证函数值是否接近零
    assert abs(f(root_document)) < 1e-2

def test_example2_newton():
    """
    验证例题2：用牛顿法求 f(x)=x^2-2 的根，x₀=1，精度 10^-4
    文档中的答案：x ≈ 1.4142
    """
    f = lambda x: x**2 - 2
    df = lambda x: 2*x
    
    # 牛顿法实现
    def newton_method(f, df, x0, tol=1e-10, max_iter=100):
        x = x0
        for i in range(max_iter):
            fx = f(x)
            if abs(fx) < tol:
                return x
            
            dfx = df(x)
            if dfx == 0:
                raise ValueError("导数为零，牛顿法失效")
            
            x_new = x - fx / dfx
            if abs(x_new - x) < tol:
                return x_new
            
            x = x_new
        
        return x
    
    # 用自己实现的牛顿法计算
    root_newton = newton_method(f, df, 1, tol=1e-4)
    # 文档中给出的答案
    root_document = 1.4142
    # 精确值是根号2
    root_exact = np.sqrt(2)
    
    # 验证自己实现的方法是否与文档答案一致
    assert is_close(root_newton, root_document, rtol=1e-4)
    # 验证文档答案是否足够接近精确值
    assert is_close(root_document, root_exact, rtol=1e-4)
    # 验证函数值是否接近零
    assert abs(f(root_document)) < 1e-4

def test_example3_newton():
    """
    验证例题3：用牛顿法求解 f(x) = e^x - x - 2 的根，初始值 x₀ = 1
    文档中的答案：x ≈ 1.146
    """
    f = lambda x: np.exp(x) - x - 2
    df = lambda x: np.exp(x) - 1
    
    # 用scipy的方法计算精确值
    root_exact = optimize.root_scalar(f, x0=1, fprime=df, method='newton').root
    # 文档中给出的答案
    root_document = 1.146
    
    # 验证文档答案是否足够接近精确值
    assert is_close(root_document, root_exact, rtol=1e-3)
    # 验证函数值是否接近零
    assert abs(f(root_document)) < 1e-2

def test_sqrt3_newton():
    """
    验证牛顿法求 sqrt(3) 的例题
    文档中的答案：sqrt(3) ≈ 1.73205
    """
    f = lambda x: x**2 - 3
    df = lambda x: 2*x
    
    # 牛顿法迭代公式的简化版本
    def newton_sqrt3(x0, iterations=3):
        x = x0
        for _ in range(iterations):
            x = 0.5 * (x + 3/x)
        return x
    
    # 用简化的牛顿法计算（文档中的迭代过程）
    root_newton = newton_sqrt3(2, iterations=3)
    # 文档中给出的答案
    root_document = 1.73205
    # 精确值是根号3
    root_exact = np.sqrt(3)
    
    # 验证自己实现的方法是否与文档答案一致
    assert is_close(root_newton, root_document, rtol=1e-5)
    # 验证文档答案是否足够接近精确值
    assert is_close(root_document, root_exact, rtol=1e-5)

def test_secant_example():
    """
    验证割线法求 f(x)=x^3+x-1 的根的例题
    文档中的答案：x ≈ 0.6822
    """
    f = lambda x: x**3 + x - 1
    
    # 割线法实现
    def secant_method(f, x0, x1, tol=1e-10, max_iter=100):
        for i in range(max_iter):
            f0, f1 = f(x0), f(x1)
            
            if abs(f1) < tol:
                return x1
            
            if f1 == f0:
                raise ValueError("割线法分母为零")
            
            x_new = x1 - f1 * (x1 - x0) / (f1 - f0)
            
            if abs(x_new - x1) < tol:
                return x_new
            
            x0, x1 = x1, x_new
        
        return x1
    
    # 用自己实现的割线法计算
    root_secant = secant_method(f, 0, 1, tol=1e-2)
    # 用scipy的方法计算精确值
    root_exact = optimize.root_scalar(f, bracket=[0, 1], method='brentq').root
    # 文档中给出的答案
    root_document = 0.6822
    
    # 验证自己实现的方法是否与文档答案一致
    assert is_close(root_secant, root_document, rtol=1e-3)
    # 验证文档答案是否足够接近精确值
    assert is_close(root_document, root_exact, rtol=1e-3)
    # 验证函数值是否接近零
    assert abs(f(root_document)) < 1e-2

def test_cos_x_example():
    """
    验证割线法求 f(x) = cos(x) - x 的根的例题
    文档中的答案：x ≈ 0.739
    """
    f = lambda x: np.cos(x) - x
    
    # 用scipy的方法计算精确值
    root_exact = optimize.root_scalar(f, bracket=[0, 1], method='brentq').root
    # 文档中给出的答案
    root_document = 0.739
    
    # 验证文档答案是否足够接近精确值
    assert is_close(root_document, root_exact, rtol=1e-3)
    # 验证函数值是否接近零
    assert abs(f(root_document)) < 1e-2

if __name__ == "__main__":
    # 运行所有测试
    pytest.main(["-v", __file__]) 