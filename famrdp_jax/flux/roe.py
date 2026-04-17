from __future__ import annotations
import jax.numpy as jnp


def _physical_flux(U, n, gamma):
    rho = jnp.maximum(U[0], 1e-300)
    u = U[1] / rho; v = U[2] / rho; w = U[3] / rho
    ke = 0.5 * (u*u + v*v + w*w)
    p = (gamma - 1.0) * (U[4] - rho * ke)
    un = u*n[0] + v*n[1] + w*n[2]
    H = (U[4] + p) / rho
    return jnp.stack([
        rho * un,
        rho * u * un + p * n[0],
        rho * v * un + p * n[1],
        rho * w * un + p * n[2],
        rho * H * un,
    ], axis=0)


def roe_flux(UL, UR, nxyz, *, gamma, eps_entropy=0.1):
    """通用法向 Roe 通量 + Harten 熵修复。"""
    rhoL = jnp.maximum(UL[0], 1e-300)
    rhoR = jnp.maximum(UR[0], 1e-300)
    uL, vL, wL = UL[1]/rhoL, UL[2]/rhoL, UL[3]/rhoL
    uR, vR, wR = UR[1]/rhoR, UR[2]/rhoR, UR[3]/rhoR
    keL = 0.5*(uL*uL+vL*vL+wL*wL); keR = 0.5*(uR*uR+vR*vR+wR*wR)
    pL = (gamma-1.0)*(UL[4]-rhoL*keL); pR = (gamma-1.0)*(UR[4]-rhoR*keR)
    HL = (UL[4]+pL)/rhoL; HR = (UR[4]+pR)/rhoR

    # Roe 平均
    sqL = jnp.sqrt(rhoL); sqR = jnp.sqrt(rhoR)
    ss = sqL + sqR
    u_h = (sqL*uL + sqR*uR)/ss
    v_h = (sqL*vL + sqR*vR)/ss
    w_h = (sqL*wL + sqR*wR)/ss
    H_h = (sqL*HL + sqR*HR)/ss
    q2  = u_h**2 + v_h**2 + w_h**2
    c2  = (gamma-1.0)*(H_h - 0.5*q2)
    c_h = jnp.sqrt(jnp.maximum(c2, 1e-30))
    un_h = u_h*nxyz[0] + v_h*nxyz[1] + w_h*nxyz[2]
    rho_h = sqL * sqR

    # Harten 修复
    eps = eps_entropy * (jnp.abs(un_h) + c_h)
    def harten(lam):
        return jnp.where(jnp.abs(lam) >= eps,
                         jnp.abs(lam),
                         0.5*(lam**2/eps + eps))

    lam1 = un_h - c_h; lam2 = un_h; lam5 = un_h + c_h
    al1 = harten(lam1); al2 = jnp.abs(lam2); al5 = harten(lam5)

    # 特征波强度 (通用法向 Roe 分解,Toro 2009 §11.3)
    dU = UR - UL
    dp = (gamma-1.0)*(dU[4] - u_h*dU[1] - v_h*dU[2] - w_h*dU[3] + 0.5*q2*dU[0])
    dun = (dU[1]*nxyz[0] + dU[2]*nxyz[1] + dU[3]*nxyz[2]) - un_h*dU[0]
    a1 = (dp - rho_h*c_h*dun) / (2.0*c_h**2)
    a5 = (dp + rho_h*c_h*dun) / (2.0*c_h**2)
    a23 = dU[0] - dp/c_h**2  # 接触 + 剪切

    # 右特征向量
    R1 = jnp.stack([jnp.ones_like(u_h), u_h-c_h*nxyz[0], v_h-c_h*nxyz[1], w_h-c_h*nxyz[2], H_h-c_h*un_h])
    R5 = jnp.stack([jnp.ones_like(u_h), u_h+c_h*nxyz[0], v_h+c_h*nxyz[1], w_h+c_h*nxyz[2], H_h+c_h*un_h])
    R2 = jnp.stack([jnp.ones_like(u_h), u_h, v_h, w_h, 0.5*q2])  # 接触波

    # 切向剪切波 (法向分量已包含在 a23 R2,需减去)
    # 完整展开:剪切速度扰动 = dU的切向动量变化
    du_t = dU[1]/rho_h - u_h*dU[0]/rho_h  # du = d(rho u)/rho - u d(rho)/rho
    dv_t = dU[2]/rho_h - v_h*dU[0]/rho_h
    dw_t = dU[3]/rho_h - w_h*dU[0]/rho_h
    # 切向速度扰动 = 全速度扰动 - 法向分量
    dun_full = du_t*nxyz[0] + dv_t*nxyz[1] + dw_t*nxyz[2]
    dvt = jnp.stack([jnp.zeros_like(u_h),
                     rho_h*(du_t - dun_full*nxyz[0]),
                     rho_h*(dv_t - dun_full*nxyz[1]),
                     rho_h*(dw_t - dun_full*nxyz[2]),
                     rho_h*(u_h*(du_t-dun_full*nxyz[0]) + v_h*(dv_t-dun_full*nxyz[1]) + w_h*(dw_t-dun_full*nxyz[2]))
                    ])

    absA_dU = al1*a1*R1 + al5*a5*R5 + al2*(a23*R2 + dvt)

    FL = _physical_flux(UL, nxyz, gamma)
    FR = _physical_flux(UR, nxyz, gamma)
    return 0.5*(FL + FR) - 0.5*absA_dU
