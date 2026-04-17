"""生成 2 块 (4,3,2) 和 (3,3,2) 的网格 fixture。"""
import numpy as np
from pathlib import Path

OUT = Path(__file__).parent / "mini.grd"

def main():
    blocks = [(4, 3, 2), (3, 3, 2)]
    with OUT.open("wb") as f:
        np.array([len(blocks)], dtype="<i4").tofile(f)
        for (ni, nj, nk) in blocks:
            np.array([ni, nj, nk], dtype="<i4").tofile(f)
        for (ni, nj, nk) in blocks:
            xs = np.linspace(0.0, 1.0, ni)
            ys = np.linspace(0.0, 1.0, nj)
            zs = np.linspace(0.0, 1.0, nk)
            X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
            for arr in (X, Y, Z):
                arr.astype("<f8").flatten(order="F").tofile(f)

if __name__ == "__main__":
    main()
    print(f"wrote {OUT}")
