# euler1d_jax: 全局 1D Euler 方程 JAX 求解器
# 数值方法：MUSCL-2 重构 + Roe 通量 + TVD-RK3
# 边界条件：Wall（滑移）/ Symmetry / Farfield（简化）
# 数据布局：所有块展平为单一全局 1D 数组 prime[N_total, 5]
