from __future__ import annotations
import numpy as np


def compare(jax_arr, ref_arr, *, name, atol=0.0, rtol=0.0, topn=10):
    assert jax_arr.shape == ref_arr.shape, f"{name}: shape {jax_arr.shape} vs {ref_arr.shape}"
    diff = np.abs(jax_arr - ref_arr)
    denom = np.maximum(np.abs(ref_arr), 1e-300)
    rel = diff / denom
    flat_idx = np.argsort(diff, axis=None)[::-1][:topn]
    coords = [np.unravel_index(i, diff.shape) for i in flat_idx]
    top = [(c, float(jax_arr[c]), float(ref_arr[c]), float(diff[c])) for c in coords]
    return {
        "ok": (diff.max() <= atol) and (rtol == 0.0 or rel.max() <= rtol),
        "max_abs": float(diff.max()),
        "max_rel": float(rel.max()),
        "top": top,
        "name": name,
    }


def report(res):
    lines = [f"[{res['name']}] max_abs={res['max_abs']:.3e} max_rel={res['max_rel']:.3e} ok={res['ok']}"]
    for (c, a, b, d) in res["top"][:5]:
        lines.append(f"  idx={c}  jax={a:.12e}  ref={b:.12e}  diff={d:.3e}")
    return "\n".join(lines)
