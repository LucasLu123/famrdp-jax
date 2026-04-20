"""度量系数计算：将 3D 坐标数组转换为 1D 全局布局的度量数据。

输出的 1D 度量数组：
    jac_flat  : (N_total,)     Jacobian（行列式）
    kxyz_flat : (N_total, 3, 3)  kxyz[l, m] = ∂ξ_l/∂x_m
    vol_flat  : (N_total,)     单元体积 = |jac|

对于 ghost 层，度量值用最近真实层的值填充（与坐标 padding 保持一致）。
"""
from __future__ import annotations
import numpy as np
import jax.numpy as jnp
from famrdp_jax.mesh.metric import compute_metrics
from .domain import DomainData, BlockInfo
from famrdp_jax.mesh.io_grd import read_grd


def build_metrics_flat(domain: DomainData, grd_path) -> tuple:
    """计算所有块的度量，展平为全局 1D 数组。

    Returns
    -------
    jac_flat  : jnp.ndarray (N_total,)
    kxyz_flat : jnp.ndarray (N_total, 3, 3)
    vol_flat  : jnp.ndarray (N_total,)
    """
    xyz_blocks = read_grd(grd_path)
    n_total = domain.n_total
    ghost   = domain.ghost

    jac_flat  = np.zeros(n_total, dtype=np.float64)
    kxyz_flat = np.zeros((n_total, 3, 3), dtype=np.float64)
    vol_flat  = np.zeros(n_total, dtype=np.float64)

    for bi in domain.blocks:
        nb = bi.block_id
        xyz_b = jnp.array(xyz_blocks[nb])   # (3, ni, nj, nk)
        g = ghost
        ni, nj, nk = bi.ni, bi.nj, bi.nk

        # 将真实格坐标 pad 到含 ghost 的缓冲区再计算度量
        buf = jnp.zeros((3, bi.ni_g, bi.nj_g, bi.nk_g), dtype=jnp.float64)
        buf = buf.at[:, g:g+ni, g:g+nj, g:g+nk].set(xyz_b)

        # edge padding: 每个轴依次向两端广播
        for ax in [1, 2, 3]:
            sl_lo = [slice(None)] * 4
            sl_hi = [slice(None)] * 4
            sl_lo[ax] = slice(0, g)
            sl_hi[ax] = slice(buf.shape[ax] - g, None)
            ref_lo = [slice(None)] * 4
            ref_hi = [slice(None)] * 4
            ref_lo[ax] = g
            ref_hi[ax] = buf.shape[ax] - g - 1
            buf = buf.at[tuple(sl_lo)].set(
                jnp.expand_dims(buf[tuple(ref_lo)], axis=ax)
            )
            buf = buf.at[tuple(sl_hi)].set(
                jnp.expand_dims(buf[tuple(ref_hi)], axis=ax)
            )

        # 计算度量（4 阶中心差，与 famrdp_jax 对齐）
        met = compute_metrics(buf, ghost=g)

        # 展平：(ni_g, nj_g, nk_g) → (N_block,)
        n_b = bi.n_total
        jac_b  = np.array(met.jac).reshape(-1)     # (N_block,)
        kxyz_b = np.array(met.kxyz).transpose(2, 3, 4, 0, 1).reshape(n_b, 3, 3)
        vol_b  = np.array(met.vol).reshape(-1)

        sl = slice(bi.offset, bi.offset + n_b)
        jac_flat[sl]    = jac_b
        kxyz_flat[sl]   = kxyz_b
        vol_flat[sl]    = vol_b

    # ghost 层度量用最近真实层填充，避免 NaN/零导致面法向量无效
    for bi in domain.blocks:
        g  = bi.ghost
        ni, nj, nk = bi.ni, bi.nj, bi.nk
        n_b = bi.n_total
        sl = slice(bi.offset, bi.offset + n_b)

        kx = kxyz_flat[sl].reshape(bi.ni_g, bi.nj_g, bi.nk_g, 3, 3)
        vl = vol_flat[sl].reshape(bi.ni_g, bi.nj_g, bi.nk_g)

        # i 方向
        kx[:g]      = kx[g:g+1];      vl[:g]      = vl[g:g+1]
        kx[g+ni:]   = kx[g+ni-1:g+ni]; vl[g+ni:]  = vl[g+ni-1:g+ni]
        # j 方向
        kx[:, :g]     = kx[:, g:g+1];     vl[:, :g]     = vl[:, g:g+1]
        kx[:, g+nj:]  = kx[:, g+nj-1:g+nj]; vl[:, g+nj:] = vl[:, g+nj-1:g+nj]
        # k 方向
        kx[:, :, :g]    = kx[:, :, g:g+1];    vl[:, :, :g]    = vl[:, :, g:g+1]
        kx[:, :, g+nk:] = kx[:, :, g+nk-1:g+nk]; vl[:, :, g+nk:] = vl[:, :, g+nk-1:g+nk]

        # 写回（reshape 返回 view 时原地修改；显式赋值保证安全）
        kxyz_flat[sl] = kx.reshape(n_b, 3, 3)
        vol_flat[sl]  = vl.reshape(n_b)

    return jnp.array(jac_flat), jnp.array(kxyz_flat), jnp.array(vol_flat)
