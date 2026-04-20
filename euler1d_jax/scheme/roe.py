"""Roe 通量求解器：(N_faces, 5) 批量接口版本。

直接移植自 famrdp_jax/flux/roe.py，适配全局 1D 的
(N_faces, 5) 原始变量接口（而非 (5, ...) 接口）。

函数签名：
    roe_flux_1d(pL, pR, normals, gamma, eps_entropy) -> F
    pL, pR  : (N_faces, 5)  原始变量 [rho, u, v, w, p]
    normals : (N_faces, 3)  面法向量（含面积权重）
    F       : (N_faces, 5)  数值通量
"""
from __future__ import annotations
import jax.numpy as jnp


def _prim_to_cons_batch(p: jnp.ndarray, gamma: float) -> jnp.ndarray:
    """(N, 5) 原始变量 → 守恒变量。"""
    rho = p[:, 0]; u = p[:, 1]; v = p[:, 2]; w = p[:, 3]; pr = p[:, 4]
    ke  = 0.5 * (u*u + v*v + w*w)
    rhoE = pr / (gamma - 1.0) + rho * ke
    return jnp.stack([rho, rho*u, rho*v, rho*w, rhoE], axis=1)


def _physical_flux_batch(U: jnp.ndarray, n: jnp.ndarray, gamma: float) -> jnp.ndarray:
    """(N, 5) 守恒变量 + (N, 3) 法向 → (N, 5) 物理法向通量。"""
    rho  = jnp.maximum(U[:, 0], 1e-300)
    u    = U[:, 1] / rho; v = U[:, 2] / rho; w = U[:, 3] / rho
    ke   = 0.5 * (u*u + v*v + w*w)
    p    = (gamma - 1.0) * (U[:, 4] - rho * ke)
    un   = u * n[:, 0] + v * n[:, 1] + w * n[:, 2]
    H    = (U[:, 4] + p) / rho
    return jnp.stack([
        rho * un,
        rho * u * un + p * n[:, 0],
        rho * v * un + p * n[:, 1],
        rho * w * un + p * n[:, 2],
        rho * H * un,
    ], axis=1)   # (N, 5)


def roe_flux_1d(
    pL: jnp.ndarray,
    pR: jnp.ndarray,
    normals: jnp.ndarray,
    *,
    gamma: float = 1.4,
    eps_entropy: float = 0.1,
) -> jnp.ndarray:
    """批量法向 Roe 通量 + Harten 熵修复。

    Parameters
    ----------
    pL, pR  : (N_faces, 5) 原始变量 [rho, u, v, w, p]
    normals : (N_faces, 3) 面法向量（已含面积，无需归一化）
    gamma   : 比热比
    eps_entropy : Harten 修复参数

    Returns
    -------
    flux : (N_faces, 5) 数值通量
    """
    # 法向归一化（仅用于 Roe 矩阵，物理通量用原始含面积法向）
    n_mag = jnp.linalg.norm(normals, axis=1, keepdims=True)   # (N, 1)
    n_hat = normals / jnp.maximum(n_mag, 1e-300)              # (N, 3) 单位法向

    # 原始量提取
    rhoL = jnp.maximum(pL[:, 0], 1e-300)
    rhoR = jnp.maximum(pR[:, 0], 1e-300)
    uL = pL[:, 1]; vL = pL[:, 2]; wL = pL[:, 3]; prL = pL[:, 4]
    uR = pR[:, 1]; vR = pR[:, 2]; wR = pR[:, 3]; prR = pR[:, 4]

    keL  = 0.5 * (uL*uL + vL*vL + wL*wL)
    keR  = 0.5 * (uR*uR + vR*vR + wR*wR)
    HL   = (prL / (gamma - 1.0) + rhoL * keL + prL) / rhoL   # H = (rhoE + p) / rho
    HR   = (prR / (gamma - 1.0) + rhoR * keR + prR) / rhoR

    # Roe 平均
    sqL   = jnp.sqrt(rhoL); sqR = jnp.sqrt(rhoR)
    ss    = jnp.maximum(sqL + sqR, 1e-150)
    u_h   = (sqL*uL + sqR*uR) / ss
    v_h   = (sqL*vL + sqR*vR) / ss
    w_h   = (sqL*wL + sqR*wR) / ss
    H_h   = (sqL*HL + sqR*HR) / ss
    q2    = u_h**2 + v_h**2 + w_h**2
    c2s   = jnp.maximum((gamma - 1.0) * (H_h - 0.5 * q2), 1e-30)
    c_h   = jnp.sqrt(c2s)
    un_h  = u_h*n_hat[:, 0] + v_h*n_hat[:, 1] + w_h*n_hat[:, 2]
    rho_h = sqL * sqR

    # Harten 熵修复
    eps_h = eps_entropy * (jnp.abs(un_h) + c_h)
    def harten(lam):
        return jnp.where(jnp.abs(lam) >= eps_h,
                         jnp.abs(lam),
                         0.5 * (lam**2 / eps_h + eps_h))

    lam1 = un_h - c_h; lam5 = un_h + c_h
    al1 = harten(lam1); al2 = jnp.abs(un_h); al5 = harten(lam5)

    # 守恒量差
    UL = _prim_to_cons_batch(pL, gamma)
    UR = _prim_to_cons_batch(pR, gamma)
    dU = UR - UL   # (N, 5)

    # 特征波强度
    dp  = (gamma - 1.0) * (dU[:, 4]
           - u_h*dU[:, 1] - v_h*dU[:, 2] - w_h*dU[:, 3]
           + 0.5*q2*dU[:, 0])
    dun = (dU[:, 1]*n_hat[:, 0] + dU[:, 2]*n_hat[:, 1] + dU[:, 3]*n_hat[:, 2]
           - un_h * dU[:, 0])
    a1  = (dp - rho_h * c_h * dun) / (2.0 * c2s)
    a5  = (dp + rho_h * c_h * dun) / (2.0 * c2s)
    a23 = dU[:, 0] - dp / c2s

    # 右特征向量 (N, 5)
    ones = jnp.ones_like(u_h)
    R1 = jnp.stack([ones, u_h - c_h*n_hat[:,0], v_h - c_h*n_hat[:,1],
                    w_h - c_h*n_hat[:,2], H_h - c_h*un_h], axis=1)
    R5 = jnp.stack([ones, u_h + c_h*n_hat[:,0], v_h + c_h*n_hat[:,1],
                    w_h + c_h*n_hat[:,2], H_h + c_h*un_h], axis=1)
    R2 = jnp.stack([ones, u_h, v_h, w_h, 0.5*q2], axis=1)

    # 切向剪切扰动
    du_t = dU[:, 1]/rho_h - u_h*dU[:, 0]/rho_h
    dv_t = dU[:, 2]/rho_h - v_h*dU[:, 0]/rho_h
    dw_t = dU[:, 3]/rho_h - w_h*dU[:, 0]/rho_h
    dun_full = du_t*n_hat[:,0] + dv_t*n_hat[:,1] + dw_t*n_hat[:,2]
    zeros = jnp.zeros_like(u_h)
    dvt = jnp.stack([
        zeros,
        rho_h*(du_t - dun_full*n_hat[:,0]),
        rho_h*(dv_t - dun_full*n_hat[:,1]),
        rho_h*(dw_t - dun_full*n_hat[:,2]),
        rho_h*(u_h*(du_t-dun_full*n_hat[:,0]) +
               v_h*(dv_t-dun_full*n_hat[:,1]) +
               w_h*(dw_t-dun_full*n_hat[:,2])),
    ], axis=1)

    # |A| ΔU = Σ αᵢ |λᵢ| Rᵢ
    absA_dU = (al1[:, None] * a1[:, None] * R1
             + al5[:, None] * a5[:, None] * R5
             + al2[:, None] * (a23[:, None] * R2 + dvt))

    # 物理通量（用含面积法向，不用归一化法向）
    FL = _physical_flux_batch(UL, normals, gamma)
    FR = _physical_flux_batch(UR, normals, gamma)

    # Roe 通量：0.5*(FL+FR) - 0.5*|A|ΔU * n_mag（矩阵乘以面积）
    # 注意：FL/FR 已用含面积法向计算，|A|ΔU 用单位法向，需乘回面积
    return 0.5 * (FL + FR) - 0.5 * absA_dU * n_mag
