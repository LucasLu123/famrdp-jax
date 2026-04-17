# xyz Reference Data — needs_reference

This directory should contain Fortran-generated reference binary files:
  block_0001.bin, block_0002.bin, ..., block_NNNN.bin

## How to generate

1. Compile the Fortran solver with `-DPROBE_XYZ` flag:
   ```
   cmake . -B build-probe -DPROBE=ON -DPROBE_WHICH=XYZ
   cmake --build build-probe
   ```
2. Create the probe output directory:
   ```
   mkdir -p probes/xyz
   ```
3. Run the solver on test3 (just the grid-read phase is sufficient):
   ```
   cd build-probe && mpirun -n 1 ./famrdp ../param.inp
   ```
4. Copy probe output here:
   ```
   cp probes/xyz/block_*.bin validation/references/test3/xyz/
   ```

## File format

Each `block_XXXX.bin` is a Fortran unformatted stream file:
- int32 x 3: ni, nj, nk
- float64 x (ni*nj*nk): x coordinates (Fortran column-major)
- float64 x (ni*nj*nk): y coordinates (Fortran column-major)
- float64 x (ni*nj*nk): z coordinates (Fortran column-major)
