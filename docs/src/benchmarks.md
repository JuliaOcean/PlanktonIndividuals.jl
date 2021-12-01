# [Benchmarks](@id benchmarks)

Here we benchmark the model performance in two `Architecture`s.
The number of individuals used in the benchmark are `(2^5, 2^10, 2^15, 2^17)`.
And we also use different grid resolutions in 2-Dimensional and 3-Dimensional model setup.

## 0-Dimensional model

This is a benchmark of a simple 0-Dimensional model setup without advection of Eulerian tracers. However, the advection of individuals still take the same amount of time whether the velocity field is provided or not.

```julia
PlanktonIndividuals v0.4.2
Julia Version 1.7.0-rc1
Commit 9eade6195e (2021-09-12 06:45 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-12.0.1 (ORCJIT, broadwell)
  GPU: Tesla P100-PCIE-12GB
```

| Arch |      N |        min |     median |       mean |        max |     memory | allocs |
|------|--------|------------|------------|------------|------------|------------|--------|
|  CPU |     32 | 978.736 μs |   1.062 ms |   1.114 ms |   1.745 ms | 639.39 KiB |   3377 |
|  CPU |   1024 |   3.217 ms |   3.319 ms |   3.357 ms |   4.003 ms | 639.39 KiB |   3377 |
|  CPU |  32768 |  73.551 ms |  73.612 ms |  73.955 ms |  77.018 ms | 638.91 KiB |   3346 |
|  CPU | 131072 | 297.726 ms | 298.756 ms | 300.489 ms | 316.688 ms | 638.91 KiB |   3346 |
|  GPU |     32 |   7.498 ms |   7.566 ms |   7.636 ms |   8.331 ms |   2.27 MiB |  16453 |
|  GPU |   1024 |   7.599 ms |   7.691 ms |   7.755 ms |   8.487 ms |   2.26 MiB |  16443 |
|  GPU |  32768 |   8.171 ms |   8.362 ms |   8.470 ms |   9.745 ms |   2.26 MiB |  16443 |
|  GPU | 131072 |   9.698 ms |  10.456 ms |  10.637 ms |  12.999 ms |   2.26 MiB |  16438 |

## 2-Dimensional model

This is the benchmark of a 2-Dimensional model setup with `(Ns, 1, Ns)` grid cells. Here `Ns = [32, 64, 128]`.

```julia
PlanktonIndividuals v0.4.2
Julia Version 1.7.0-rc1
Commit 9eade6195e (2021-09-12 06:45 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-12.0.1 (ORCJIT, broadwell)
  GPU: Tesla P100-PCIE-12GB
```

| Arch |      N |  Ns |        min |     median |       mean |        max |    memory | allocs |
|------|--------|-----|------------|------------|------------|------------|-----------|--------|
|  CPU |     32 |  32 |   4.183 ms |   5.037 ms |   4.932 ms |   5.125 ms |  2.86 MiB |   3413 |
|  CPU |     32 |  64 |  12.474 ms |  12.583 ms |  12.697 ms |  13.670 ms |  8.84 MiB |   3386 |
|  CPU |     32 | 128 |  46.953 ms |  57.432 ms |  53.965 ms |  61.537 ms | 31.87 MiB |   3386 |
|  CPU |   1024 |  32 |   6.800 ms |   7.802 ms |   7.589 ms |   7.937 ms |  2.86 MiB |   3413 |
|  CPU |   1024 |  64 |  15.106 ms |  15.227 ms |  15.361 ms |  16.476 ms |  8.84 MiB |   3386 |
|  CPU |   1024 | 128 |  51.023 ms |  61.336 ms |  57.659 ms |  62.330 ms | 31.87 MiB |   3386 |
|  CPU |  32768 |  32 |  91.757 ms |  91.996 ms |  92.255 ms |  93.695 ms |  2.86 MiB |   3382 |
|  CPU |  32768 |  64 | 105.509 ms | 105.603 ms | 106.028 ms | 108.820 ms |  8.84 MiB |   3386 |
|  CPU |  32768 | 128 | 154.187 ms | 155.702 ms | 156.419 ms | 163.824 ms | 31.87 MiB |   3386 |
|  CPU | 131072 |  32 | 362.675 ms | 363.038 ms | 363.071 ms | 363.607 ms |  2.86 MiB |   3382 |
|  CPU | 131072 |  64 | 392.255 ms | 392.962 ms | 395.636 ms | 405.071 ms |  8.84 MiB |   3386 |
|  CPU | 131072 | 128 | 447.502 ms | 458.867 ms | 461.654 ms | 488.007 ms | 31.87 MiB |   3386 |
|  GPU |     32 |  32 |   8.094 ms |   8.161 ms |   8.285 ms |   9.522 ms |  2.29 MiB |  16137 |
|  GPU |     32 |  64 |   7.603 ms |   7.783 ms |   7.833 ms |   8.644 ms |  2.39 MiB |  16141 |
|  GPU |     32 | 128 |   7.728 ms |   7.783 ms |   7.966 ms |   9.569 ms |  2.76 MiB |  16221 |
|  GPU |   1024 |  32 |   8.248 ms |   8.310 ms |   8.432 ms |   9.660 ms |  2.29 MiB |  16127 |
|  GPU |   1024 |  64 |   7.253 ms |   7.329 ms |   7.428 ms |   8.332 ms |  2.38 MiB |  16131 |
|  GPU |   1024 | 128 |   7.957 ms |   7.991 ms |   8.173 ms |   9.711 ms |  2.76 MiB |  16211 |
|  GPU |  32768 |  32 |   8.173 ms |   8.251 ms |   8.372 ms |   9.494 ms |  2.29 MiB |  16127 |
|  GPU |  32768 |  64 |   7.237 ms |   7.291 ms |   7.435 ms |   8.777 ms |  2.38 MiB |  16131 |
|  GPU |  32768 | 128 |   7.681 ms |   7.816 ms |   8.036 ms |  10.264 ms |  2.76 MiB |  16211 |
|  GPU | 131072 |  32 |   8.970 ms |   9.371 ms |   9.390 ms |   9.851 ms |  2.29 MiB |  16153 |
|  GPU | 131072 |  64 |   9.451 ms |  10.731 ms |  10.602 ms |  10.960 ms |  2.38 MiB |  16126 |
|  GPU | 131072 | 128 |   9.267 ms |  12.095 ms |  11.808 ms |  12.248 ms |  2.76 MiB |  16206 |

## 3-Dimensional model

This is the benchmark of a 3-Dimensional model setup with `(Ns, Ns, Ns)` grid cells. Here `Ns = [32, 64]`.

```julia
PlanktonIndividuals v0.4.2
Julia Version 1.7.0-rc1
Commit 9eade6195e (2021-09-12 06:45 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-12.0.1 (ORCJIT, broadwell)
  GPU: Tesla P100-PCIE-12GB
```

| Arch |      N |  Ns |        min |     median |       mean |        max |    memory | allocs |
|------|--------|-----|------------|------------|------------|------------|-----------|--------|
|  CPU |     32 |  32 |  38.263 ms |  38.316 ms |  39.038 ms |  41.863 ms |  1.54 MiB |   3154 |
|  CPU |     32 |  64 | 332.699 ms | 333.257 ms | 333.191 ms | 333.711 ms |  8.59 MiB |   3155 |
|  CPU |   1024 |  32 |  41.214 ms |  41.334 ms |  41.623 ms |  44.407 ms |  1.54 MiB |   3154 |
|  CPU |   1024 |  64 | 337.645 ms | 341.374 ms | 350.123 ms | 375.033 ms |  8.59 MiB |   3155 |
|  CPU |  32768 |  32 | 135.441 ms | 135.510 ms | 135.875 ms | 137.648 ms |  1.54 MiB |   3154 |
|  CPU |  32768 |  64 | 447.552 ms | 448.844 ms | 458.740 ms | 499.685 ms |  8.59 MiB |   3155 |
|  CPU | 131072 |  32 | 433.618 ms | 433.704 ms | 433.846 ms | 434.720 ms |  1.54 MiB |   3154 |
|  CPU | 131072 |  64 | 763.314 ms | 763.408 ms | 777.291 ms | 848.858 ms |  8.59 MiB |   3155 |
|  GPU |     32 |  32 |   7.094 ms |   7.159 ms |   7.348 ms |   9.046 ms |  3.26 MiB |  15561 |
|  GPU |     32 |  64 |  10.841 ms |  11.494 ms |  11.443 ms |  11.617 ms | 10.31 MiB |  15611 |
|  GPU |   1024 |  32 |   6.679 ms |   6.790 ms |   6.897 ms |   8.001 ms |  3.25 MiB |  15551 |
|  GPU |   1024 |  64 |  10.791 ms |  11.485 ms |  11.427 ms |  11.617 ms | 10.30 MiB |  15601 |
|  GPU |  32768 |  32 |   6.686 ms |   6.762 ms |   6.936 ms |   8.584 ms |  3.25 MiB |  15551 |
|  GPU |  32768 |  64 |  11.470 ms |  11.857 ms |  11.821 ms |  12.028 ms | 10.30 MiB |  15601 |
|  GPU | 131072 |  32 |   8.724 ms |  10.342 ms |  10.180 ms |  10.585 ms |  3.25 MiB |  15546 |
|  GPU | 131072 |  64 |  12.760 ms |  15.537 ms |  15.228 ms |  15.779 ms | 10.30 MiB |  15627 |
