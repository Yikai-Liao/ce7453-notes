# CE7453 章节练习题验证

这个目录包含用于验证CE7453课程各章节练习题答案正确性的Python脚本。

## 环境要求

- Python 3.6+
- NumPy
- SciPy
- pytest

## 安装依赖

```bash
pip install numpy scipy pytest
```

## 使用方法

1. 进入`verification`目录
   ```bash
   cd verification
   ```

2. 运行所有测试
   ```bash
   pytest
   ```

3. 运行特定章节的测试
   ```bash
   pytest test_chapter01.py
   ```

## 测试结果

| 章节 | 测试结果 | 备注 |
|------|---------|------|
| 第1章 - 根查找方法 | ✅ 全部通过 | 二分法、牛顿法和割线法的例题答案都正确 |
| 第2章 - 线性系统 | ✅ 全部通过 | 高斯消元法、LU分解、迭代法例题答案都正确 |
| 第3章 - Bezier/B-spline/插值 | ✅ 全部通过 | Bezier曲线、B-spline曲线和插值例题答案都正确 |

## 文件结构

- `test_chapter01.py`: 第1章 - 根查找方法的验证
- `test_chapter02.py`: 第2章 - 线性系统的验证
- `test_chapter03.py`: 第3章 - Bezier和B样条插值的验证
- `test_chapter04.py`: 第4章 - 数值微分和积分的验证（待完成）
- `test_chapter05.py`: 第5章 - 最小二乘法的验证（待完成）
- `test_chapter06.py`: 第6章 - 特征值分析的验证（待完成）
- `test_chapter07.py`: 第7章 - 傅立叶变换的验证（待完成）
- `test_chapter08.py`: 第8章 - 3D数据配准的验证（待完成）

## 测试说明

每个测试文件对应课程中的一个章节，包含对该章节中练习题答案的验证测试。当所有测试通过时，表示章节中的练习题答案是正确的。

## 修复和注意事项

- 第2章中的Gauss-Seidel迭代法例题计算过程在原文中是正确的
- 第3章中的B-spline曲线计算也需要按照原文中的方式逐步计算 