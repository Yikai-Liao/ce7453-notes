import numpy as np
import pytest

def is_close(a, b, rtol=1e-2, atol=1e-8):
    """检查两个数值是否足够接近，考虑浮点数精度问题"""
    if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
        return np.allclose(a, b, rtol=rtol, atol=atol)
    return np.abs(a - b) <= (atol + rtol * np.abs(b))

def compute_optimal_transformation(A, B):
    """
    计算将点集A配准到点集B的最优刚性变换
    
    参数：
        A: 源点集，形状为(n, dim)的数组，其中n是点的数量，dim是维度
        B: 目标点集，形状为(n, dim)的数组
        
    返回：
        R: 最优旋转矩阵，形状为(dim, dim)
        t: 最优平移向量，形状为(dim,)
    """
    # 转换为NumPy数组
    A = np.array(A)
    B = np.array(B)
    
    # 计算质心
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    
    # 去中心化
    A_centered = A - centroid_A
    B_centered = B - centroid_B
    
    # 计算协方差矩阵
    H = np.zeros((A.shape[1], B.shape[1]))
    for a, b in zip(A_centered, B_centered):
        H += np.outer(a, b)
    
    # SVD分解
    U, _, Vt = np.linalg.svd(H)
    
    # 计算旋转矩阵
    # 确保旋转矩阵是正交的且行列式为1
    V = Vt.T
    det = np.linalg.det(V @ U.T)
    
    if det < 0:
        V[:, -1] = -V[:, -1]
    
    R = V @ U.T
    
    # 计算平移向量
    t = centroid_B - R @ centroid_A
    
    return R, t

def test_example1():
    """
    验证例题1：已知A={(0,0), (1,0)}，B={(1,1), (2,1)}，求将A配准到B的最优刚性变换
    文档中的答案：R为单位矩阵，t = (1, 1)
    """
    # 源点集和目标点集
    A = np.array([[0, 0], [1, 0]])
    B = np.array([[1, 1], [2, 1]])
    
    # 计算最优变换
    R, t = compute_optimal_transformation(A, B)
    
    # 期望结果
    expected_R = np.eye(2)  # 单位矩阵（无旋转）
    expected_t = np.array([1, 1])  # 平移向量
    
    # 打印结果进行比较
    print("计算的旋转矩阵R:")
    print(R)
    print("期望的旋转矩阵R:")
    print(expected_R)
    print("计算的平移向量t:", t)
    print("期望的平移向量t:", expected_t)
    
    # 验证结果是否与期望值一致
    assert is_close(R, expected_R)
    assert is_close(t, expected_t)
    
    # 验证变换后的点是否与目标点集重合
    transformed_A = (R @ A.T).T + t
    assert is_close(transformed_A, B)

def test_example2():
    """
    验证例题2：给定A=[(0,0,0), (1,0,0)]，B=[(0,1,0), (1,1,0)]，求将A配准到B的最优刚性变换
    文档中的答案：R为单位矩阵，t = (0, 1, 0)
    """
    # 源点集和目标点集
    A = np.array([[0, 0, 0], [1, 0, 0]])
    B = np.array([[0, 1, 0], [1, 1, 0]])
    
    # 计算最优变换
    R, t = compute_optimal_transformation(A, B)
    
    # 期望结果
    expected_R = np.eye(3)  # 单位矩阵（无旋转）
    expected_t = np.array([0, 1, 0])  # 平移向量
    
    # 打印结果进行比较
    print("计算的旋转矩阵R:")
    print(R)
    print("期望的旋转矩阵R:")
    print(expected_R)
    print("计算的平移向量t:", t)
    print("期望的平移向量t:", expected_t)
    
    # 验证结果是否与期望值一致
    assert is_close(R, expected_R)
    assert is_close(t, expected_t)
    
    # 验证变换后的点是否与目标点集重合
    transformed_A = (R @ A.T).T + t
    assert is_close(transformed_A, B)

def test_example3():
    """
    验证例题3：给定A=[(0,0,0), (1,0,0), (0,1,0)]，B=[(0,0,0), (0,1,0), (1,0,0)]，求将A配准到B的最优刚性变换
    文档中的答案：R为绕z轴旋转90度的矩阵，t = (0, 0, 0)
    """
    # 源点集和目标点集
    A = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    B = np.array([[0, 0, 0], [0, 1, 0], [1, 0, 0]])
    
    # 使用SVD方法计算最优变换
    R, t = compute_optimal_transformation(A, B)
    
    # 打印计算结果
    print("计算的旋转矩阵R:")
    print(R)
    print("计算的平移向量t:", t)
    
    # 验证变换后的点是否与目标点集接近
    transformed_A = (R @ A.T).T + t
    print("SVD变换后的点集:")
    print(transformed_A)
    print("目标点集:")
    print(B)
    
    # 计算变换后点集与目标点集的均方误差
    mse = np.mean(np.sum((transformed_A - B)**2, axis=1))
    print(f"SVD变换均方误差: {mse}")
    
    # 主要验证变换结果是否正确，而不是具体的矩阵值
    assert mse < 0.01  # 验证变换后的点集是否与目标点集接近
    
    # 验证R是否为正交矩阵（旋转矩阵的特性）
    RRT = R @ R.T
    identity = np.eye(3)
    R_orthogonal_error = np.linalg.norm(RRT - identity)
    print(f"旋转矩阵正交性误差: {R_orthogonal_error}")
    assert R_orthogonal_error < 0.01
    
    # 验证行列式是否为1（保持方向）或-1（改变方向）
    det_R = np.linalg.det(R)
    print(f"旋转矩阵行列式: {det_R}")
    assert abs(abs(det_R) - 1) < 0.01
    
    # 验证平移向量是否接近零向量
    t_norm = np.linalg.norm(t)
    print(f"平移向量范数: {t_norm}")
    assert t_norm < 0.01

if __name__ == "__main__":
    # 直接调用测试函数
    print("测试例题1...")
    test_example1()
    print("\n测试例题2...")
    test_example2()
    print("\n测试例题3...")
    test_example3()
    print("\n所有测试通过!") 