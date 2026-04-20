"""TVD-RK3 时间积分（Shu-Osher 三步法）。

仅对真实单元更新；ghost 层在每个 RK 子步前由 BC/halo 重新设置。

守恒变量更新策略：
    - 时间推进在原始变量 prime 上进行（直接更新 prime）
    - 对于 Euler 方程，RHS 本身是原始变量的函数，TVD-RK3 可直接作用于原始变量
    - 物理约束：rho > 0, p > 0，必要时裁剪

Shu-Osher 系数：
    U1 = U0 + dt * L(U0)
    U2 = 3/4 * U0 + 1/4 * (U1 + dt * L(U1))
    U_new = 1/3 * U0 + 2/3 * (U2 + dt * L(U2))
"""
from __future__ import annotations
from functools import partial
import jax
import jax.numpy as jnp

from euler1d_jax.bc.ghost import apply_bc_all, apply_halo
from euler1d_jax.solver.rhs import compute_euler_rhs


def step_euler(
    prime: jnp.ndarray,
    precomp: dict,
    bc_ops: list,
    cut_maps: list,
    true_idx: jnp.ndarray,
    dt: float,
    *,
    gamma: float = 1.4,
    limiter: str = "minmod",
) -> jnp.ndarray:
    """推进一步（TVD-RK3）。

    Parameters
    ----------
    prime    : (N_total, 5) 当前原始变量（含 ghost 层）
    precomp  : 预计算索引/法向量
    bc_ops   : BC 操作列表（build_bc_ops 输出）
    cut_maps : CUT1TO1 映射列表（build_cut1to1_maps 输出）
    true_idx : (N_true,) 真实单元全局索引
    dt       : 时间步长
    gamma    : 比热比
    limiter  : MUSCL-2 限制器

    Returns
    -------
    prime_new : (N_total, 5) 推进后原始变量
    """
    def L(p):
        """RHS 算子：先施加 BC，再计算 Euler RHS。"""
        p = apply_halo(p, cut_maps)
        p = apply_bc_all(p, bc_ops)
        return compute_euler_rhs(p, precomp, gamma=gamma, limiter=limiter)

    # Stage 1
    rhs0 = L(prime)
    prime1 = prime.at[true_idx].add(dt * rhs0[true_idx])
    prime1 = _clip_physical(prime1, true_idx)

    # Stage 2
    rhs1 = L(prime1)
    prime2_true = (0.75 * prime[true_idx]
                   + 0.25 * (prime1[true_idx] + dt * rhs1[true_idx]))
    prime2 = prime.at[true_idx].set(prime2_true)
    prime2 = _clip_physical(prime2, true_idx)

    # Stage 3
    rhs2 = L(prime2)
    prime_new_true = ((1.0 / 3.0) * prime[true_idx]
                      + (2.0 / 3.0) * (prime2[true_idx] + dt * rhs2[true_idx]))
    prime_new = prime.at[true_idx].set(prime_new_true)
    prime_new = _clip_physical(prime_new, true_idx)

    return prime_new


def make_step_jit(
    precomp: dict,
    bc_ops: list,
    cut_maps: list,
    true_idx,
    *,
    gamma: float = 1.4,
    limiter: str = "minmod",
):
    """返回 JIT 编译的单步推进函数 step(prime, dt) -> prime_new。

    所有预计算数据（precomp/bc_ops/cut_maps/true_idx）通过闭包捕获为常量，
    JAX 编译时展开 Python 循环和分支，生成融合的 XLA 计算图。
    首次调用时触发编译（约数秒），后续调用直接执行 XLA 代码。
    """
    @jax.jit
    def _step(prime: jnp.ndarray, dt: float) -> jnp.ndarray:
        return step_euler(prime, precomp, bc_ops, cut_maps, true_idx, dt,
                          gamma=gamma, limiter=limiter)
    return _step


def _clip_physical(prime: jnp.ndarray, true_idx: jnp.ndarray) -> jnp.ndarray:
    """裁剪真实单元的 rho 和 p，防止非物理负值。"""
    rho_min = 1e-10
    p_min   = 1e-10
    p_true = prime[true_idx]
    p_true = p_true.at[:, 0].set(jnp.maximum(p_true[:, 0], rho_min))
    p_true = p_true.at[:, 4].set(jnp.maximum(p_true[:, 4], p_min))
    return prime.at[true_idx].set(p_true)


def compute_dt_cfl(
    prime: jnp.ndarray,
    vol_flat: jnp.ndarray,
    kxyz_flat: jnp.ndarray,
    true_idx: jnp.ndarray,
    cfl: float,
    *,
    gamma: float = 1.4,
) -> float:
    """基于 CFL 条件计算全局最小时间步长。

    谱半径：r_d = |u_n| + c  (d 方向法向速度 + 声速)
    dt = CFL * min(vol / (r_0 + r_1 + r_2))
    """
    p = prime[true_idx]          # (N_true, 5)
    k = kxyz_flat[true_idx]     # (N_true, 3, 3)
    v = vol_flat[true_idx]       # (N_true,)

    rho = jnp.maximum(p[:, 0], 1e-300)
    u   = p[:, 1]; vel = p[:, 2]; w = p[:, 3]; pr = p[:, 4]
    uvw = jnp.stack([u, vel, w], axis=1)    # (N_true, 3)
    c   = jnp.sqrt(gamma * pr / rho)        # (N_true,)

    sr_total = jnp.zeros(p.shape[0], dtype=jnp.float64)
    for d in range(3):
        # k[:, d, :] = ∂ξ_d/∂x，法向 = kxyz[d, :]
        n_d   = k[:, d, :]                  # (N_true, 3) 协变基向量
        n_mag = jnp.linalg.norm(n_d, axis=1)        # (N_true,) 面积
        n_hat = n_d / jnp.maximum(n_mag, 1e-300)[:, None]    # 单位法向
        vn    = jnp.sum(uvw * n_hat, axis=1)        # 法向速度
        sr_total += (jnp.abs(vn) + c) * n_mag

    dt_local = cfl * v / jnp.maximum(sr_total, 1e-300)
    return float(jnp.min(dt_local))
