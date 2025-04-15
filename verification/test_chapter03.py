import numpy as np
import pytest
from scipy import interpolate

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

def test_bezier_de_casteljau():
    """
    验证例题：已知控制点 P₀=(0,0), P₁=(1,2), P₂=(3,3)，求 t=0.5 时的 Bezier 曲线点
    文档中的答案：P(0.5) = (1.25, 1.75)
    """
    # 控制点
    P0 = np.array([0, 0])
    P1 = np.array([1, 2])
    P2 = np.array([3, 3])
    t = 0.5
    
    # 直接手动计算，按照文档中的步骤
    P0_1 = (1 - t) * P0 + t * P1
    P1_1 = (1 - t) * P1 + t * P2
    P0_2 = (1 - t) * P0_1 + t * P1_1
    
    # 验证手动计算结果
    # P0_1 = (1-0.5)(0,0) + 0.5(1,2) = (0.5,1)
    assert is_array_close(P0_1, np.array([0.5, 1.0]))
    # P1_1 = (1-0.5)(1,2) + 0.5(3,3) = (2,2.5)
    assert is_array_close(P1_1, np.array([2.0, 2.5]))
    # P0_2 = (1-0.5)(0.5,1) + 0.5(2,2.5) = (1.25,1.75)
    assert is_array_close(P0_2, np.array([1.25, 1.75]))
    
    # 文档中给出的答案
    expected_point = np.array([1.25, 1.75])
    
    # 验证最终结果是否与文档答案一致
    assert is_array_close(P0_2, expected_point)

def test_bezier_additional():
    """
    验证新增例题：计算给定控制点 P₀=(0,0)，P₁=(1,1)，P₂=(2,0) 的二次Bezier曲线在 t=0.5 处的点
    文档中的答案：(1, 0.5)
    """
    # 控制点
    P0 = np.array([0, 0])
    P1 = np.array([1, 1])
    P2 = np.array([2, 0])
    t = 0.5
    
    # 手动计算
    P0_1 = (1 - t) * P0 + t * P1
    P1_1 = (1 - t) * P1 + t * P2
    P0_2 = (1 - t) * P0_1 + t * P1_1
    
    # 文档中给出的答案
    expected_point = np.array([1, 0.5])
    
    # 验证手动计算是否与文档答案一致
    assert is_array_close(P0_2, expected_point)

def test_bspline():
    """
    验证例题：已知三次 B-spline 曲线，控制点 P₀=(0,0), P₁=(1,2), P₂=(3,3), P₃=(4,0)，
    节点向量 [0,0,0,0,1,2,2,2,2]，求 u=1 处的曲线点
    文档中的答案：r(1) = (2.75, 2.0)
    """
    # 控制点
    control_points = np.array([[0, 0], [1, 2], [3, 3], [4, 0]])
    
    # 按照文档中的计算，直接计算基函数值和最终点
    # 文档中的基函数值
    N14 = 0.25  # N_{1,4}(1) = 0.25
    N24 = 0.5   # N_{2,4}(1) = 0.5
    N34 = 0.25  # N_{3,4}(1) = 0.25
    
    # 手动计算曲线点
    manual_point = N14 * control_points[1] + N24 * control_points[2] + N34 * control_points[3]
    
    # 文档中给出的答案
    expected_point = np.array([2.75, 2.0])
    
    # 验证手动计算是否与文档答案一致
    assert is_array_close(manual_point, expected_point)
    
    # 检查手动计算过程
    P1 = control_points[1]  # (1,2)
    P2 = control_points[2]  # (3,3)
    P3 = control_points[3]  # (4,0)
    
    # r(1) = 0.25*(1,2) + 0.5*(3,3) + 0.25*(4,0)
    #      = (0.25*1 + 0.5*3 + 0.25*4, 0.25*2 + 0.5*3 + 0.25*0)
    #      = (0.25 + 1.5 + 1, 0.5 + 1.5 + 0)
    #      = (2.75, 2.0)
    
    # 验证计算步骤
    x_value = 0.25*1 + 0.5*3 + 0.25*4
    y_value = 0.25*2 + 0.5*3 + 0.25*0
    
    assert is_close(x_value, 2.75)
    assert is_close(y_value, 2.0)
    assert is_array_close(np.array([x_value, y_value]), expected_point)

def test_lagrange_interpolation():
    """
    验证例题：已知点 (0,1), (1,2), (2,0)，用 Lagrange 插值求 x=1.5 处的函数值
    文档中的答案：L(1.5) = 1.375
    """
    # 数据点
    x = np.array([0, 1, 2])
    y = np.array([1, 2, 0])
    
    # 插值点
    x_interp = 1.5
    
    # Lagrange 插值实现
    def lagrange_interpolation(x_data, y_data, x_interp):
        n = len(x_data)
        result = 0.0
        for i in range(n):
            li = 1.0
            for j in range(n):
                if i != j:
                    li *= (x_interp - x_data[j]) / (x_data[i] - x_data[j])
            result += y_data[i] * li
        return result
    
    # 使用Lagrange插值计算
    interp_value = lagrange_interpolation(x, y, x_interp)
    
    # 文档中给出的答案
    expected_value = 1.375
    
    # 验证计算结果是否与文档答案一致
    assert is_close(interp_value, expected_value)
    
    # 使用scipy的Lagrange多项式验证
    lagrange_poly = interpolate.lagrange(x, y)
    scipy_value = lagrange_poly(x_interp)
    
    # 验证scipy计算结果是否与文档答案一致
    assert is_close(scipy_value, expected_value)

def test_newton_interpolation():
    """
    验证例题：已知点 (0,1), (1,2), (2,0)，用 Newton 插值求 x=1.5 处的函数值
    文档中的答案：N(1.5) = 1.375
    """
    # 数据点
    x = np.array([0, 1, 2])
    y = np.array([1, 2, 0])
    
    # 插值点
    x_interp = 1.5
    
    # 计算差商表
    def divided_differences(x, y):
        n = len(x)
        coef = np.zeros([n, n])
        # 填充第一列，即函数值
        coef[:, 0] = y
        
        for j in range(1, n):
            for i in range(n - j):
                coef[i, j] = (coef[i + 1, j - 1] - coef[i, j - 1]) / (x[i + j] - x[i])
                
        return coef[0, :]  # 返回第一行作为系数
    
    # Newton 插值实现
    def newton_interpolation(x_data, y_data, x_interp):
        coef = divided_differences(x_data, y_data)
        n = len(x_data)
        result = coef[0]
        
        for i in range(1, n):
            term = coef[i]
            for j in range(i):
                term *= (x_interp - x_data[j])
            result += term
            
        return result
    
    # 使用Newton插值计算
    interp_value = newton_interpolation(x, y, x_interp)
    
    # 文档中给出的答案
    expected_value = 1.375
    
    # 验证计算结果是否与文档答案一致
    assert is_close(interp_value, expected_value)
    
    # 也可以检查使用其他方法计算的插值结果
    poly = np.polyfit(x, y, 2)  # 拟合二次多项式
    polyval = np.polyval(poly, x_interp)
    
    # 验证多项式拟合结果是否与预期一致
    assert is_close(polyval, expected_value)

if __name__ == "__main__":
    # 运行所有测试
    pytest.main(["-v", __file__]) 