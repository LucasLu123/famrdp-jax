# validation/fortran_probes/

Fortran probe 代码,用于把内部量 dump 为二进制供 JAX 对拍。

## 启用

```bash
cmake . -B build-probe -DPROBE=ON -DPROBE_WHICH=METRIC,RHS,STEP
cmake --build build-probe -j
```

## 输出

- probes/metric/block_NNNN.bin
- probes/rhs/block_NNNN.bin  
- probes/step/step_NNNNNN.bin

## 格式

Fortran unformatted stream,含 shape 头 + 数据,Python 用 numpy.fromfile 读。
