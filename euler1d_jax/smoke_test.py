"""快速 smoke test: 验证预处理和一步推进正确运行。"""
import sys
import pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
try:
    jax.config.update("jax_default_device", jax.devices("gpu")[0])
except RuntimeError:
    jax.config.update("jax_default_device", jax.devices("cpu")[0])
import jax.numpy as jnp

from euler1d_jax.mesh.domain import build_domain, build_cut1to1_maps
from euler1d_jax.mesh.metrics import build_metrics_flat
from euler1d_jax.scheme.precomp import build_precomp
from euler1d_jax.bc.ghost import build_bc_ops, apply_bc_all, apply_halo
from euler1d_jax.solver.rhs import compute_euler_rhs
from euler1d_jax.solver.rk3 import step_euler, compute_dt_cfl

GRD = "1d_famrdp/grid_files/Cylinder/test3.grd"
INP = "1d_famrdp/grid_files/Cylinder/test3.inp"
Q_INF = np.array([1.0, 0.1, 0.0, 0.0, 0.7142857142857143])
GAMMA = 1.4

print("1. 构建域索引...", flush=True)
domain   = build_domain(GRD, INP, ghost=2)
cut_maps = build_cut1to1_maps(domain)
true_idx = jnp.array(domain.true_idx)
print(f"   blocks={len(domain.blocks)}, n_total={domain.n_total}, true={len(domain.true_idx)}")

print("2. 计算度量系数...", flush=True)
jac, kxyz, vol = build_metrics_flat(domain, GRD)
print(f"   vol: min={float(vol.min()):.4e}, max={float(vol.max()):.4e}")
assert float(vol[true_idx].min()) > 0, "真实单元体积应为正"

print("3. 生成面索引...", flush=True)
precomp = build_precomp(domain, kxyz, vol)
for d in range(3):
    print(f"   dir={d}: {precomp[d]['stencil'].shape[1]} faces")

print("4. 构建 BC 操作表...", flush=True)
bc_ops = build_bc_ops(domain, Q_INF)
print(f"   {len(bc_ops)} BC operations")

print("5. 初始化流场...", flush=True)
prime = jnp.broadcast_to(jnp.array(Q_INF)[None, :], (domain.n_total, 5)).copy()
prime = apply_halo(prime, cut_maps)
prime = apply_bc_all(prime, bc_ops)
print(f"   prime: shape={prime.shape}, dtype={prime.dtype}")

print("6. 计算单步 RHS...", flush=True)
rhs = compute_euler_rhs(prime, precomp, gamma=GAMMA)
print(f"   rhs: max={float(jnp.max(jnp.abs(rhs[true_idx]))):.4e}")
assert jnp.isfinite(rhs).all(), "RHS 含 NaN/inf"

print("7. 计算 CFL 时间步...", flush=True)
dt = compute_dt_cfl(prime, vol, kxyz, true_idx, cfl=0.5, gamma=GAMMA)
print(f"   dt = {dt:.4e}")
assert dt > 0 and dt < 1.0, f"dt 不合理: {dt}"

print("8. 推进 5 步...", flush=True)
import time
t0 = time.perf_counter()
for i in range(5):
    prime_prev = prime
    prime = step_euler(prime, precomp, bc_ops, cut_maps, true_idx, dt,
                       gamma=GAMMA, limiter="minmod")
    jax.block_until_ready(prime)
    res = float(jnp.max(jnp.abs(prime[true_idx] - prime_prev[true_idx])))
    print(f"   step {i+1}: residual={res:.4e}", flush=True)
t1 = time.perf_counter()
print(f"   5步耗时: {t1-t0:.2f}s, avg={((t1-t0)/5)*1000:.1f}ms/step")

assert jnp.isfinite(prime).all(), "最终流场含 NaN/inf"
print("\n=== SMOKE TEST PASSED ===")
