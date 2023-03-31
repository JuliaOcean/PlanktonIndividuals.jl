# [Benchmarks](@id benchmarks)

Here we benchmark the model performance in two `Architecture`s.
The number of individuals used in the benchmark are `(2^10, 2^15, 2^17, 2^20)`.
And we also use different grid resolutions in 2-Dimensional and 3-Dimensional model setup.

## 0-Dimensional model

This is a benchmark of a simple 0-Dimensional model setup without advection of Eulerian tracers. However, the advection of individuals still take the same amount of time whether the velocity field is provided or not.

```julia
PlanktonIndividuals v0.6.1
Julia Version 1.8.0
Commit 5544a0fab76 (2022-08-17 13:38 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-13.0.1 (ORCJIT, broadwell)
  GPU: Tesla P100-PCIE-12GB
  CUDA runtime 11.8, artifact installation
  CUDA driver 11.2
  NVIDIA driver 460.84.0
```

| Arch |       N |        min |     median |       mean |        max |     memory | allocs |
|------|---------|------------|------------|------------|------------|------------|--------|
|  CPU |    1024 |   2.945 ms |   3.016 ms |   3.167 ms |   4.328 ms | 478.67 KiB |   2992 |
|  CPU |   32768 |  69.741 ms |  69.812 ms |  71.594 ms |  80.231 ms | 477.72 KiB |   2931 |
|  CPU |  131072 | 276.553 ms | 276.966 ms | 280.569 ms | 300.907 ms | 477.72 KiB |   2931 |
|  CPU | 1048576 |    2.582 s |    2.590 s |    2.590 s |    2.598 s | 477.72 KiB |   2931 |
|  GPU |    1024 |   7.085 ms |   7.158 ms |   7.364 ms |   9.323 ms |   1.92 MiB |  21327 |
|  GPU |   32768 |   7.435 ms |   7.520 ms |   7.925 ms |  10.173 ms |   1.92 MiB |  21327 |
|  GPU |  131072 |   7.053 ms |   9.161 ms |   9.851 ms |  19.812 ms |   1.92 MiB |  21294 |
|  GPU | 1048576 |   8.005 ms |  46.217 ms |  47.484 ms | 122.516 ms |   1.92 MiB |  21294 |

## 2-Dimensional model

This is the benchmark of a 2-Dimensional model setup with `(Ns, 1, Ns)` grid cells. Here `Ns = [32, 64, 128]`.

```julia
PlanktonIndividuals v0.6.1
Julia Version 1.8.0
Commit 5544a0fab76 (2022-08-17 13:38 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-13.0.1 (ORCJIT, broadwell)
  GPU: Tesla P100-PCIE-12GB
  CUDA runtime 11.8, artifact installation
  CUDA driver 11.2
  NVIDIA driver 460.84.0
```

| Arch |       N |  Ns |        min |     median |       mean |        max |    memory | allocs |
|------|---------|-----|------------|------------|------------|------------|-----------|--------|
|  CPU |    1024 |  32 |   8.096 ms |   8.132 ms |   8.211 ms |   8.688 ms |  2.70 MiB |   3109 |
|  CPU |    1024 |  64 |  19.889 ms |  19.940 ms |  20.064 ms |  20.952 ms |  8.68 MiB |   3052 |
|  CPU |    1024 | 128 |  68.735 ms |  69.030 ms |  69.672 ms |  75.046 ms | 31.72 MiB |   3052 |
|  CPU |   32768 |  32 |  74.115 ms |  74.154 ms |  76.313 ms |  85.288 ms |  2.70 MiB |   3048 |
|  CPU |   32768 |  64 |  89.999 ms |  90.163 ms |  92.340 ms | 101.475 ms |  8.68 MiB |   3052 |
|  CPU |   32768 | 128 | 162.286 ms | 162.618 ms | 168.129 ms | 190.011 ms | 31.72 MiB |   3052 |
|  CPU |  131072 |  32 | 282.810 ms | 282.913 ms | 286.631 ms | 307.620 ms |  2.70 MiB |   3048 |
|  CPU |  131072 |  64 | 328.584 ms | 328.962 ms | 332.448 ms | 357.787 ms |  8.68 MiB |   3052 |
|  CPU |  131072 | 128 | 447.271 ms | 453.263 ms | 470.108 ms | 509.040 ms | 31.72 MiB |   3052 |
|  CPU | 1048576 |  32 |    2.476 s |    2.476 s |    2.501 s |    2.552 s |  2.70 MiB |   3048 |
|  CPU | 1048576 |  64 |    2.910 s |    2.911 s |    2.911 s |    2.911 s |  8.68 MiB |   3052 |
|  CPU | 1048576 | 128 |    2.905 s |    2.909 s |    2.909 s |    2.914 s | 31.72 MiB |   3052 |
|  GPU |    1024 |  32 |   6.902 ms |   6.920 ms |   7.101 ms |   8.719 ms |  1.98 MiB |  21513 |
|  GPU |    1024 |  64 |   7.417 ms |   7.622 ms |   7.755 ms |   8.430 ms |  2.07 MiB |  21632 |
|  GPU |    1024 | 128 |   7.734 ms |   8.071 ms |   8.141 ms |   8.854 ms |  2.45 MiB |  21713 |
|  GPU |   32768 |  32 |   7.011 ms |   7.092 ms |   7.392 ms |  10.142 ms |  1.98 MiB |  21513 |
|  GPU |   32768 |  64 |   6.769 ms |   6.837 ms |   7.152 ms |  10.035 ms |  2.07 MiB |  21632 |
|  GPU |   32768 | 128 |   7.027 ms |   8.381 ms |   8.561 ms |  11.845 ms |  2.45 MiB |  21713 |
|  GPU |  131072 |  32 |   6.580 ms |   8.054 ms |   8.560 ms |  15.323 ms |  1.98 MiB |  21541 |
|  GPU |  131072 |  64 |   7.491 ms |   9.106 ms |   9.664 ms |  16.128 ms |  2.07 MiB |  21599 |
|  GPU |  131072 | 128 |   7.918 ms |  12.640 ms |  12.791 ms |  23.534 ms |  2.45 MiB |  21680 |
|  GPU | 1048576 |  32 |   9.781 ms |  35.539 ms |  36.437 ms |  59.171 ms |  1.98 MiB |  21528 |
|  GPU | 1048576 |  64 |  10.682 ms |  37.958 ms |  39.055 ms |  65.476 ms |  2.08 MiB |  21647 |
|  GPU | 1048576 | 128 |   7.994 ms |  50.094 ms |  50.772 ms | 126.537 ms |  2.45 MiB |  21680 |

## 3-Dimensional model

This is the benchmark of a 3-Dimensional model setup with `(Ns, Ns, Ns)` grid cells. Here `Ns = [32, 64]`.

```julia
PlanktonIndividuals v0.6.1
Julia Version 1.8.0
Commit 5544a0fab76 (2022-08-17 13:38 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-13.0.1 (ORCJIT, broadwell)
  GPU: Tesla P100-PCIE-12GB
  CUDA runtime 11.8, artifact installation
  CUDA driver 11.2
  NVIDIA driver 460.84.0
```

| Arch |       N |  Ns |        min |     median |       mean |        max |   memory | allocs |
|------|---------|-----|------------|------------|------------|------------|----------|--------|
|  CPU |    1024 |  32 |  50.081 ms |  50.249 ms |  50.421 ms |  51.994 ms | 1.38 MiB |   2820 |
|  CPU |    1024 |  64 | 410.840 ms | 459.105 ms | 451.043 ms | 459.516 ms | 8.43 MiB |   2821 |
|  CPU |   32768 |  32 | 124.176 ms | 124.312 ms | 126.438 ms | 138.224 ms | 1.38 MiB |   2820 |
|  CPU |   32768 |  64 | 498.713 ms | 534.237 ms | 534.148 ms | 554.501 ms | 8.43 MiB |   2821 |
|  CPU |  131072 |  32 | 351.282 ms | 351.674 ms | 355.733 ms | 387.071 ms | 1.38 MiB |   2820 |
|  CPU |  131072 |  64 | 790.994 ms | 808.337 ms | 816.691 ms | 848.149 ms | 8.43 MiB |   2821 |
|  CPU | 1048576 |  32 |    3.019 s |    3.072 s |    3.072 s |    3.125 s | 1.38 MiB |   2820 |
|  CPU | 1048576 |  64 |    3.258 s |    3.258 s |    3.258 s |    3.258 s | 8.43 MiB |   2821 |
|  GPU |    1024 |  32 |   6.229 ms |   6.286 ms |   6.466 ms |   7.329 ms | 2.94 MiB |  21053 |
|  GPU |    1024 |  64 |   9.194 ms |  11.891 ms |  11.689 ms |  12.604 ms | 9.99 MiB |  21077 |
|  GPU |   32768 |  32 |   6.570 ms |   6.638 ms |   6.966 ms |   8.974 ms | 2.94 MiB |  21053 |
|  GPU |   32768 |  64 |   9.143 ms |  12.882 ms |  12.712 ms |  15.781 ms | 9.99 MiB |  21077 |
|  GPU |  131072 |  32 |   6.481 ms |   9.150 ms |   9.469 ms |  16.907 ms | 2.94 MiB |  21081 |
|  GPU |  131072 |  64 |   9.212 ms |  16.623 ms |  16.438 ms |  25.557 ms | 9.99 MiB |  21105 |
|  GPU | 1048576 |  32 |   7.257 ms |  39.894 ms |  40.268 ms |  96.189 ms | 2.94 MiB |  21020 |
|  GPU | 1048576 |  64 |   9.586 ms |  54.934 ms |  53.741 ms | 118.675 ms | 9.99 MiB |  21105 |
