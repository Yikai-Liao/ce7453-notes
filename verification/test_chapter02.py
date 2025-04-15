import numpy as np
import pytest
from scipy import linalg
import sys
import os

# 直接定义需要的函数
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

def check_equation_solution(A, x, b, tol=1e-8):
    """
    检查线性方程Ax=b的解是否正确
    """
    residual = A @ x - b
    return np.linalg.norm(residual) < tol

def test_example1_gaussian():
    """
    验证例题1：用高斯消元法解线性方程组
    2x + 3y = 8
    5x + 4y = 13
    文档中的答案：x = 1, y = 2
    """
    A = np.array([[2, 3], [5, 4]])
    b = np.array([8, 13])
    
    # 用numpy的solve函数计算精确解
    x_exact = np.linalg.solve(A, b)
    # 文档中给出的答案
    x_document = np.array([1, 2])
    
    # 验证文档答案是否足够接近精确解
    assert is_array_close(x_document, x_exact)
    # 验证Ax=b是否满足
    assert check_equation_solution(A, x_document, b)

def test_example2_gaussian():
    """
    验证例题2：用高斯消元法解线性方程组
    3x + 2y = 7
    x + 4y = 5
    文档中的答案：x = 1.8, y = 0.8
    """
    A = np.array([[3, 2], [1, 4]])
    b = np.array([7, 5])
    
    # 用numpy的solve函数计算精确解
    x_exact = np.linalg.solve(A, b)
    # 文档中给出的答案
    x_document = np.array([1.8, 0.8])
    
    # 验证文档答案是否足够接近精确解
    assert is_array_close(x_document, x_exact)
    # 验证Ax=b是否满足
    assert check_equation_solution(A, x_document, b)

def test_example3_lu():
    """
    验证例题3：用LU分解解线性方程组
    4x + 3y = 24
    8x + 7y = 52
    文档中的答案：x = 3, y = 4
    """
    A = np.array([[4, 3], [8, 7]])
    b = np.array([24, 52])
    
    # 使用LU分解求解
    def solve_lu(A, b):
        # LU分解
        P, L, U = linalg.lu(A)
        # 排列b
        b_perm = P @ b
        # 前代求解Ly = b_perm
        y = linalg.solve_triangular(L, b_perm, lower=True)
        # 回代求解Ux = y
        x = linalg.solve_triangular(U, y, lower=False)
        return x
    
    # 用LU分解方法计算解
    x_lu = solve_lu(A, b)
    # 文档中给出的答案
    x_document = np.array([3, 4])
    
    # 验证文档答案是否足够接近精确解
    assert is_array_close(x_document, x_lu)
    # 验证Ax=b是否满足
    assert check_equation_solution(A, x_document, b)
    
    # 检查文档中的L和U矩阵是否正确
    L_document = np.array([[1, 0], [2, 1]])
    U_document = np.array([[4, 3], [0, 1]])
    # 验证A = LU
    assert is_array_close(A, L_document @ U_document)

def test_example4_jacobi():
    """
    验证例题4：用Jacobi迭代法解线性方程组
    4x + y = 9
    x + 3y = 7
    初始值 x^(0) = 0, y^(0) = 0，迭代两次
    文档中的答案：x^(2) = 1.67, y^(2) = 1.58
    """
    A = np.array([[4, 1], [1, 3]])
    b = np.array([9, 7])
    
    # Jacobi迭代法实现
    def jacobi_iteration(A, b, x0, iterations):
        n = len(b)
        x = x0.copy()
        for k in range(iterations):
            x_new = np.zeros(n)
            for i in range(n):
                s = 0
                for j in range(n):
                    if j != i:
                        s += A[i, j] * x[j]
                x_new[i] = (b[i] - s) / A[i, i]
            x = x_new
        return x
    
    # 用Jacobi迭代法计算
    x0 = np.array([0, 0])
    x_jacobi = jacobi_iteration(A, b, x0, 2)
    
    # 文档中给出的答案
    x_document = np.array([1.67, 1.58])
    
    # 验证文档答案是否足够接近计算结果
    assert is_array_close(x_document, x_jacobi, rtol=1e-2)
    
    # 精确解
    x_exact = np.linalg.solve(A, b)
    # 检查迭代结果是否向精确解靠近
    assert np.linalg.norm(x_jacobi - x_exact) < np.linalg.norm(x0 - x_exact)

def test_example5_gauss_seidel():
    """
    验证例题5：用Gauss-Seidel迭代法解线性方程组
    4x + y = 9
    x + 3y = 7
    初始值 x^(0) = 0, y^(0) = 0，迭代两次
    文档中的答案：x^(2) = 1.86, y^(2) = 1.71
    """
    A = np.array([[4, 1], [1, 3]])
    b = np.array([9, 7])
    
    # 这个测试直接按照文档中的计算过程逐步计算
    # 文档的Gauss-Seidel迭代过程如下：
    
    # 初始值
    x0, y0 = 0, 0
    
    # 第一次迭代
    x1 = (9 - y0) / 4  # = 9/4 = 2.25
    y1 = (7 - x1) / 3  # = (7 - 2.25) / 3 = 1.58333...
    
    # 第二次迭代
    x2 = (9 - y1) / 4  # = (9 - 1.58333...) / 4 = 1.85416...
    y2 = (7 - x2) / 3  # = (7 - 1.85416...) / 3 = 1.71527...
    
    # 与文档中的结果比较（文档结果经过四舍五入）
    assert abs(x2 - 1.86) < 0.01
    assert abs(y2 - 1.71) < 0.01
    
    # 或者使用我们的工具函数检查
    assert is_close(x2, 1.86, rtol=1e-2)
    assert is_close(y2, 1.71, rtol=1e-2)
    
    # 确认文档中的答案
    doc_result = np.array([1.86, 1.71])
    calculated_result = np.array([x2, y2])
    assert is_array_close(doc_result, calculated_result, rtol=1e-2)
    
    # 精确解
    x_exact = np.linalg.solve(A, b)  # [1.5, 1.5]
    
    # 验证迭代是否朝着精确解收敛
    assert np.linalg.norm(calculated_result - x_exact) < np.linalg.norm(np.array([x0, y0]) - x_exact)

def test_example6_gaussian_4x4():
    """
    验证例题6：用高斯消元法解4阶线性方程组(修正版)
    2x₁ + 1x₂ + 3x₃ + 2x₄ = 21
    4x₁ + 3x₂ + 2x₃ + 1x₄ = 20
    1x₁ + 2x₂ + 4x₃ + 3x₄ = 29
    3x₁ + 4x₂ + 1x₃ + 2x₄ = 22
    文档修正后的答案：x₁ = 1, x₂ = 2, x₃ = 3, x₄ = 4
    """
    A = np.array([
        [2, 1, 3, 2],
        [4, 3, 2, 1],
        [1, 2, 4, 3],
        [3, 4, 1, 2]
    ])
    b = np.array([21, 20, 29, 22])
    
    # 用numpy的solve函数计算精确解
    x_exact = np.linalg.solve(A, b)
    # 文档中给出的答案
    x_document = np.array([1, 2, 3, 4])
    
    # 验证文档答案是否足够接近精确解
    assert is_array_close(x_document, x_exact)
    # 验证Ax=b是否满足
    assert check_equation_solution(A, x_document, b)

def test_example7_jacobi_convergence():
    """
    验证例题7：判断迭代法收敛性
    10x + 2y = 12
    3x + 15y = 18
    文档中的结论：Jacobi迭代法对于该线性方程组是收敛的。
    """
    A = np.array([[10, 2], [3, 15]])
    
    # 检查是否严格对角占优
    def is_strictly_diagonally_dominant(A):
        n = A.shape[0]
        for i in range(n):
            diagonal = abs(A[i, i])
            row_sum = sum(abs(A[i, j]) for j in range(n) if j != i)
            if diagonal <= row_sum:
                return False
        return True
    
    # 检查谱半径条件
    def check_spectral_radius(A):
        D = np.diag(np.diag(A))
        D_inv = np.linalg.inv(D)
        M = -D_inv @ (A - D)
        eigenvalues = np.linalg.eigvals(M)
        return max(abs(eigenvalues)) < 1
    
    # 验证是否严格对角占优
    assert is_strictly_diagonally_dominant(A)
    # 验证迭代矩阵的谱半径是否小于1
    assert check_spectral_radius(A)

if __name__ == "__main__":
    # 运行所有测试
    pytest.main(["-v", __file__]) 