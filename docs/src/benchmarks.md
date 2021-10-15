# [Benchmarks](@id benchmarks)

Here we benchmark the model performance in two `Architecture`s.
The number of individuals used in the benchmark are `(2, 2^5, 2^10, 2^15)`.
And we also use different grid resolutions in 2-Dimensional and 3-Dimensional model setup.

## 0-Dimensional model

This is a benchmark of a simple 0-Dimensional model setup without advection of Eulerian tracers. However, the advection of individuals still take the same amount of time whether the velocity field is provided or not.

```julia
PlanktonIndividuals v0.3.5
Julia Version 1.6.2
Commit 1b93d53fc4 (2021-07-14 15:36 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, broadwell)
  GPU: Tesla P100-PCIE-12GB
```

| Arch |     N |        min |       mean |        max |     memory | allocs |
|------|-------|------------|------------|------------|------------|--------|
|  CPU |     2 | 891.520 μs | 993.516 μs |   1.527 ms | 630.64 KiB |   3403 |
|  CPU |    32 | 966.725 μs |   1.325 ms |   3.304 ms | 630.64 KiB |   3403 |
|  CPU |  1024 |   3.187 ms |   3.332 ms |   3.999 ms | 630.64 KiB |   3403 |
|  CPU | 32768 |  73.771 ms |  76.160 ms |  81.231 ms | 630.16 KiB |   3372 |
|  GPU |     2 |   8.287 ms |   8.512 ms |   9.272 ms |   2.51 MiB |  26341 |
|  GPU |    32 |   8.386 ms |   8.520 ms |   8.888 ms |   2.51 MiB |  26321 |
|  GPU |  1024 |   8.335 ms |   8.546 ms |   9.378 ms |   2.51 MiB |  26332 |
|  GPU | 32768 |   8.484 ms |   8.645 ms |   9.274 ms |   2.51 MiB |  26332 |

## 2-Dimensional model

This is the benchmark of a 2-Dimensional model setup with `(Ns, 1, Ns)` grid cells. Here `Ns = [32, 64, 128]`.

```julia
PlanktonIndividuals v0.3.5
Julia Version 1.6.2
Commit 1b93d53fc4 (2021-07-14 15:36 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, broadwell)
  GPU: Tesla P100-PCIE-12GB
```

| Arch |     N |  Ns |        min |       mean |        max |    memory | allocs |
|------|-------|-----|------------|------------|------------|-----------|--------|
|  CPU |     2 |  32 |   4.333 ms |   4.788 ms |   5.391 ms |  2.85 MiB |   3439 |
|  CPU |     2 |  64 |  13.344 ms |  17.459 ms |  18.877 ms |  8.83 MiB |   3412 |
|  CPU |     2 | 128 |  52.593 ms |  60.354 ms |  67.061 ms | 31.87 MiB |   3412 |
|  CPU |    32 |  32 |   4.457 ms |   5.318 ms |   5.529 ms |  2.85 MiB |   3439 |
|  CPU |    32 |  64 |  13.408 ms |  15.995 ms |  17.091 ms |  8.83 MiB |   3412 |
|  CPU |    32 | 128 |  53.254 ms |  55.793 ms |  72.153 ms | 31.87 MiB |   3412 |
|  CPU |  1024 |  32 |   7.104 ms |   7.796 ms |   8.288 ms |  2.85 MiB |   3439 |
|  CPU |  1024 |  64 |  17.687 ms |  20.254 ms |  22.347 ms |  8.83 MiB |   3412 |
|  CPU |  1024 | 128 |  54.135 ms |  55.020 ms |  61.095 ms | 31.87 MiB |   3412 |
|  CPU | 32768 |  32 |  93.750 ms |  93.982 ms |  95.164 ms |  2.85 MiB |   3408 |
|  CPU | 32768 |  64 | 107.631 ms | 109.786 ms | 113.385 ms |  8.83 MiB |   3412 |
|  CPU | 32768 | 128 | 146.828 ms | 147.768 ms | 153.719 ms | 31.87 MiB |   3412 |
|  GPU |     2 |  32 |   8.190 ms |   8.363 ms |   9.151 ms |  2.56 MiB |  26026 |
|  GPU |     2 |  64 |   8.197 ms |   8.370 ms |   9.112 ms |  2.65 MiB |  26109 |
|  GPU |     2 | 128 |   8.286 ms |   8.518 ms |   9.901 ms |  3.03 MiB |  26254 |
|  GPU |    32 |  32 |   8.069 ms |   8.230 ms |   8.798 ms |  2.56 MiB |  25998 |
|  GPU |    32 |  64 |   8.240 ms |   8.425 ms |   9.406 ms |  2.65 MiB |  26101 |
|  GPU |    32 | 128 |   8.370 ms |   8.603 ms |   9.845 ms |  3.03 MiB |  26216 |
|  GPU |  1024 |  32 |   8.109 ms |   8.231 ms |   8.743 ms |  2.56 MiB |  26017 |
|  GPU |  1024 |  64 |   9.132 ms |   9.514 ms |  10.792 ms |  2.65 MiB |  26114 |
|  GPU |  1024 | 128 |   8.320 ms |   8.527 ms |   9.971 ms |  3.03 MiB |  26237 |
|  GPU | 32768 |  32 |   8.159 ms |   8.320 ms |   9.133 ms |  2.56 MiB |  26037 |
|  GPU | 32768 |  64 |   8.591 ms |   8.814 ms |   9.488 ms |  2.65 MiB |  26136 |
|  GPU | 32768 | 128 |   8.324 ms |   8.554 ms |   9.934 ms |  3.03 MiB |  26233 |

## 3-Dimensional model

This is the benchmark of a 3-Dimensional model setup with `(Ns, Ns, Ns)` grid cells. Here `Ns = [32, 64]`.

```julia
PlanktonIndividuals v0.3.5
Julia Version 1.6.2
Commit 1b93d53fc4 (2021-07-14 15:36 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, broadwell)
  GPU: Tesla P100-PCIE-12GB
```

| Arch |     N |  Ns |        min |       mean |        max |   memory  | allocs |
|------|-------|-----|------------|------------|------------|-----------|--------|
|  CPU |     2 |  32 |  38.397 ms |  39.306 ms |  45.373 ms |  1.53 MiB |   3180 |
|  CPU |     2 |  64 | 326.354 ms | 330.216 ms | 343.420 ms |  8.58 MiB |   3181 |
|  CPU |    32 |  32 |  38.239 ms |  38.611 ms |  41.134 ms |  1.53 MiB |   3180 |
|  CPU |    32 |  64 | 326.684 ms | 328.294 ms | 334.325 ms |  8.58 MiB |   3181 |
|  CPU |  1024 |  32 |  41.330 ms |  42.156 ms |  46.322 ms |  1.53 MiB |   3180 |
|  CPU |  1024 |  64 | 329.427 ms | 329.809 ms | 330.225 ms |  8.58 MiB |   3181 |
|  CPU | 32768 |  32 | 136.383 ms | 136.740 ms | 138.306 ms |  1.53 MiB |   3180 |
|  CPU | 32768 |  64 | 479.461 ms | 483.093 ms | 511.492 ms |  8.58 MiB |   3181 |
|  GPU |     2 |  32 |   7.783 ms |   7.967 ms |   9.286 ms |  3.51 MiB |  25128 |
|  GPU |     2 |  64 |  11.031 ms |  11.614 ms |  14.124 ms | 10.57 MiB |  25668 |
|  GPU |    32 |  32 |   8.651 ms |   9.927 ms |  17.521 ms |  3.51 MiB |  25108 |
|  GPU |    32 |  64 |  10.505 ms |  10.758 ms |  12.446 ms | 10.57 MiB |  25654 |
|  GPU |  1024 |  32 |   7.803 ms |   8.307 ms |   9.483 ms |  3.51 MiB |  25155 |
|  GPU |  1024 |  64 |  10.502 ms |  10.755 ms |  12.111 ms | 10.57 MiB |  25693 |
|  GPU | 32768 |  32 |   7.808 ms |   8.043 ms |   9.871 ms |  3.51 MiB |  25153 |
|  GPU | 32768 |  64 |  11.533 ms |  12.037 ms |  13.389 ms | 10.57 MiB |  25745 |
