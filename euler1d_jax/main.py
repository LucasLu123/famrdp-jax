"""euler1d_jax 主程序：1D 全局布局 Euler 方程求解器。

用法：
    cd euler1d_jax
    conda run -n jax python main.py [param.json]
"""
from __future__ import annotations
import sys
import os
import json
import time
import pathlib

# 把 famrdp_jax 加入搜索路径（与 euler1d_jax 在同一父目录）
_HERE = pathlib.Path(__file__).parent
sys.path.insert(0, str(_HERE.parent))

import numpy as np

# JAX 配置
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
from euler1d_jax.solver.rk3 import step_euler, compute_dt_cfl, make_step_jit


def load_param(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def make_uniform_init(domain, q_inf: np.ndarray) -> jnp.ndarray:
    """用自由流值初始化所有单元（含 ghost 层）。"""
    prime = np.broadcast_to(q_inf[None, :], (domain.n_total, 5)).copy()
    return jnp.array(prime)


def compute_residual(prime: jnp.ndarray, prime_prev: jnp.ndarray,
                     true_idx: jnp.ndarray) -> tuple:
    diff = jnp.abs(prime[true_idx] - prime_prev[true_idx])
    res_max   = float(jnp.max(diff))
    res_total = float(jnp.sum(diff))
    return res_total, res_max


def main():
    param_file = sys.argv[1] if len(sys.argv) > 1 else "param.json"
    p = load_param(param_file)

    # 路径解析（相对于 param.json 所在目录）
    base = pathlib.Path(param_file).parent
    grd  = base / p["grid_file"]
    inp  = base / p["top_file"]

    ghost  = p["ghost"]
    gamma  = p["gamma"]
    limiter = p["limiter"]
    inf    = p["inflow"]
    q_inf  = np.array([inf["rho"], inf["u"], inf["v"], inf["w"], inf["p"]])

    n_steps       = p["time"]["n_steps"]
    cfl           = p["time"]["cfl"]
    dtau          = p["time"]["dtau"]
    use_cfl       = p["time"]["use_cfl"]
    print_interval = p["time"]["print_interval"]

    out_dir = base / p["output"]["residual_file"].rsplit("/", 1)[0]
    out_dir.mkdir(parents=True, exist_ok=True)
    res_file = base / p["output"]["residual_file"]
    sol_file = base / p["output"]["solution_file"]

    print(f"euler1d_jax: grid={grd.name}, ghost={ghost}, limiter={limiter}")
    print(f"  Euler 方程，Ma={inf['mach']:.2f}，{n_steps} 步")

    # -----------------------------------------------------------------------
    # 1. 预处理
    # -----------------------------------------------------------------------
    t0 = time.perf_counter()
    print("正在构建域索引...", flush=True)
    domain   = build_domain(grd, inp, ghost=ghost)
    cut_maps = build_cut1to1_maps(domain)
    true_idx = jnp.array(domain.true_idx)
    print(f"  总单元数（含 ghost）：{domain.n_total}，真实单元：{len(domain.true_idx)}")

    print("正在计算度量系数...", flush=True)
    jac_flat, kxyz_flat, vol_flat = build_metrics_flat(domain, grd)

    print("正在生成面索引和法向量...", flush=True)
    precomp = build_precomp(domain, kxyz_flat, vol_flat)

    print("正在构建 BC 操作表...", flush=True)
    bc_ops = build_bc_ops(domain, q_inf)

    t1 = time.perf_counter()
    print(f"预处理完成，耗时 {t1-t0:.2f}s")

    # -----------------------------------------------------------------------
    # 2. 初始化
    # -----------------------------------------------------------------------
    prime = make_uniform_init(domain, q_inf)
    prime = apply_halo(prime, cut_maps)
    prime = apply_bc_all(prime, bc_ops)

    # -----------------------------------------------------------------------
    # 3. JIT 编译（首次调用触发，预热）
    # -----------------------------------------------------------------------
    print("正在 JIT 编译 step 函数（首次编译约需数秒）...", flush=True)
    step_jit = make_step_jit(precomp, bc_ops, cut_maps, true_idx,
                             gamma=gamma, limiter=limiter)
    # 预热：用一个假步进触发编译，不计入计时
    _dt0 = compute_dt_cfl(prime, vol_flat, kxyz_flat, true_idx, cfl, gamma=gamma)
    _prime_w = step_jit(prime, _dt0)
    jax.block_until_ready(_prime_w)
    t_jit = time.perf_counter()
    print(f"JIT 编译完成，耗时 {t_jit-t1:.2f}s")

    # -----------------------------------------------------------------------
    # 4. 主循环
    # -----------------------------------------------------------------------
    print(f"\n{'STEP':>6}  {'res_total':>12}  {'res_max':>12}  {'dt':>10}")
    residuals = []
    avg_times = []
    t_start   = time.perf_counter()

    for step in range(n_steps):
        prime_prev = prime

        if use_cfl:
            dt = compute_dt_cfl(prime, vol_flat, kxyz_flat, true_idx, cfl, gamma=gamma)
        else:
            dt = dtau

        prime = step_jit(prime, dt)
        jax.block_until_ready(prime)

        t_now = time.perf_counter()
        avg_times.append(t_now - t_start)
        t_start = t_now

        res_total, res_max = compute_residual(prime, prime_prev, true_idx)
        residuals.append((step, res_total, res_max, float(dt)))

        if step % print_interval == 0:
            print(f"{step:6d}  {res_total:12.5e}  {res_max:12.5e}  {dt:10.4e}")

    # -----------------------------------------------------------------------
    # 5. 输出
    # -----------------------------------------------------------------------
    avg_step = np.mean(avg_times[1:]) if len(avg_times) > 1 else avg_times[0]
    total_t  = sum(avg_times)
    print(f"\navg: {avg_step:.6f}s/step, total_time {total_t:.2f}s")

    # 保存残差
    with open(res_file, "w") as f:
        f.write("step,res_total,res_max,dt\n")
        for row in residuals:
            f.write(f"{row[0]},{row[1]:.8e},{row[2]:.8e},{row[3]:.8e}\n")

    # 保存最终流场（真实单元）
    prime_true = np.array(prime[true_idx])
    np.savez(sol_file, prime=prime_true, true_idx=np.array(true_idx))
    print(f"结果已写入: {sol_file}")


if __name__ == "__main__":
    main()
