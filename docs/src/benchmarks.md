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

| Arch | N       | min        | median     | mean       | max        | memory     | allocs |
| ---- | ------- | ---------- | ---------- | ---------- | ---------- | ---------- | ------ |
| CPU  | 1024    | 2.945 ms   | 3.016 ms   | 3.167 ms   | 4.328 ms   | 478.67 KiB | 2992   |
| CPU  | 32768   | 69.741 ms  | 69.812 ms  | 71.594 ms  | 80.231 ms  | 477.72 KiB | 2931   |
| CPU  | 131072  | 276.553 ms | 276.966 ms | 280.569 ms | 300.907 ms | 477.72 KiB | 2931   |
| CPU  | 1048576 | 2.582 s    | 2.590 s    | 2.590 s    | 2.598 s    | 477.72 KiB | 2931   |
| GPU  | 1024    | 7.085 ms   | 7.158 ms   | 7.364 ms   | 9.323 ms   | 1.92 MiB   | 21327  |
| GPU  | 32768   | 7.435 ms   | 7.520 ms   | 7.925 ms   | 10.173 ms  | 1.92 MiB   | 21327  |
| GPU  | 131072  | 7.053 ms   | 9.161 ms   | 9.851 ms   | 19.812 ms  | 1.92 MiB   | 21294  |
| GPU  | 1048576 | 8.005 ms   | 46.217 ms  | 47.484 ms  | 122.516 ms | 1.92 MiB   | 21294  |

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

| Arch | N       | Ns  | min        | median     | mean       | max        | memory    | allocs |
| ---- | ------- | --- | ---------- | ---------- | ---------- | ---------- | --------- | ------ |
| CPU  | 1024    | 32  | 8.096 ms   | 8.132 ms   | 8.211 ms   | 8.688 ms   | 2.70 MiB  | 3109   |
| CPU  | 1024    | 64  | 19.889 ms  | 19.940 ms  | 20.064 ms  | 20.952 ms  | 8.68 MiB  | 3052   |
| CPU  | 1024    | 128 | 68.735 ms  | 69.030 ms  | 69.672 ms  | 75.046 ms  | 31.72 MiB | 3052   |
| CPU  | 32768   | 32  | 74.115 ms  | 74.154 ms  | 76.313 ms  | 85.288 ms  | 2.70 MiB  | 3048   |
| CPU  | 32768   | 64  | 89.999 ms  | 90.163 ms  | 92.340 ms  | 101.475 ms | 8.68 MiB  | 3052   |
| CPU  | 32768   | 128 | 162.286 ms | 162.618 ms | 168.129 ms | 190.011 ms | 31.72 MiB | 3052   |
| CPU  | 131072  | 32  | 282.810 ms | 282.913 ms | 286.631 ms | 307.620 ms | 2.70 MiB  | 3048   |
| CPU  | 131072  | 64  | 328.584 ms | 328.962 ms | 332.448 ms | 357.787 ms | 8.68 MiB  | 3052   |
| CPU  | 131072  | 128 | 447.271 ms | 453.263 ms | 470.108 ms | 509.040 ms | 31.72 MiB | 3052   |
| CPU  | 1048576 | 32  | 2.476 s    | 2.476 s    | 2.501 s    | 2.552 s    | 2.70 MiB  | 3048   |
| CPU  | 1048576 | 64  | 2.910 s    | 2.911 s    | 2.911 s    | 2.911 s    | 8.68 MiB  | 3052   |
| CPU  | 1048576 | 128 | 2.905 s    | 2.909 s    | 2.909 s    | 2.914 s    | 31.72 MiB | 3052   |
| GPU  | 1024    | 32  | 6.902 ms   | 6.920 ms   | 7.101 ms   | 8.719 ms   | 1.98 MiB  | 21513  |
| GPU  | 1024    | 64  | 7.417 ms   | 7.622 ms   | 7.755 ms   | 8.430 ms   | 2.07 MiB  | 21632  |
| GPU  | 1024    | 128 | 7.734 ms   | 8.071 ms   | 8.141 ms   | 8.854 ms   | 2.45 MiB  | 21713  |
| GPU  | 32768   | 32  | 7.011 ms   | 7.092 ms   | 7.392 ms   | 10.142 ms  | 1.98 MiB  | 21513  |
| GPU  | 32768   | 64  | 6.769 ms   | 6.837 ms   | 7.152 ms   | 10.035 ms  | 2.07 MiB  | 21632  |
| GPU  | 32768   | 128 | 7.027 ms   | 8.381 ms   | 8.561 ms   | 11.845 ms  | 2.45 MiB  | 21713  |
| GPU  | 131072  | 32  | 6.580 ms   | 8.054 ms   | 8.560 ms   | 15.323 ms  | 1.98 MiB  | 21541  |
| GPU  | 131072  | 64  | 7.491 ms   | 9.106 ms   | 9.664 ms   | 16.128 ms  | 2.07 MiB  | 21599  |
| GPU  | 131072  | 128 | 7.918 ms   | 12.640 ms  | 12.791 ms  | 23.534 ms  | 2.45 MiB  | 21680  |
| GPU  | 1048576 | 32  | 9.781 ms   | 35.539 ms  | 36.437 ms  | 59.171 ms  | 1.98 MiB  | 21528  |
| GPU  | 1048576 | 64  | 10.682 ms  | 37.958 ms  | 39.055 ms  | 65.476 ms  | 2.08 MiB  | 21647  |
| GPU  | 1048576 | 128 | 7.994 ms   | 50.094 ms  | 50.772 ms  | 126.537 ms | 2.45 MiB  | 21680  |

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

| Arch | N       | Ns  | min        | median     | mean       | max        | memory   | allocs |
| ---- | ------- | --- | ---------- | ---------- | ---------- | ---------- | -------- | ------ |
| CPU  | 1024    | 32  | 50.081 ms  | 50.249 ms  | 50.421 ms  | 51.994 ms  | 1.38 MiB | 2820   |
| CPU  | 1024    | 64  | 410.840 ms | 459.105 ms | 451.043 ms | 459.516 ms | 8.43 MiB | 2821   |
| CPU  | 32768   | 32  | 124.176 ms | 124.312 ms | 126.438 ms | 138.224 ms | 1.38 MiB | 2820   |
| CPU  | 32768   | 64  | 498.713 ms | 534.237 ms | 534.148 ms | 554.501 ms | 8.43 MiB | 2821   |
| CPU  | 131072  | 32  | 351.282 ms | 351.674 ms | 355.733 ms | 387.071 ms | 1.38 MiB | 2820   |
| CPU  | 131072  | 64  | 790.994 ms | 808.337 ms | 816.691 ms | 848.149 ms | 8.43 MiB | 2821   |
| CPU  | 1048576 | 32  | 3.019 s    | 3.072 s    | 3.072 s    | 3.125 s    | 1.38 MiB | 2820   |
| CPU  | 1048576 | 64  | 3.258 s    | 3.258 s    | 3.258 s    | 3.258 s    | 8.43 MiB | 2821   |
| GPU  | 1024    | 32  | 6.229 ms   | 6.286 ms   | 6.466 ms   | 7.329 ms   | 2.94 MiB | 21053  |
| GPU  | 1024    | 64  | 9.194 ms   | 11.891 ms  | 11.689 ms  | 12.604 ms  | 9.99 MiB | 21077  |
| GPU  | 32768   | 32  | 6.570 ms   | 6.638 ms   | 6.966 ms   | 8.974 ms   | 2.94 MiB | 21053  |
| GPU  | 32768   | 64  | 9.143 ms   | 12.882 ms  | 12.712 ms  | 15.781 ms  | 9.99 MiB | 21077  |
| GPU  | 131072  | 32  | 6.481 ms   | 9.150 ms   | 9.469 ms   | 16.907 ms  | 2.94 MiB | 21081  |
| GPU  | 131072  | 64  | 9.212 ms   | 16.623 ms  | 16.438 ms  | 25.557 ms  | 9.99 MiB | 21105  |
| GPU  | 1048576 | 32  | 7.257 ms   | 39.894 ms  | 40.268 ms  | 96.189 ms  | 2.94 MiB | 21020  |
| GPU  | 1048576 | 64  | 9.586 ms   | 54.934 ms  | 53.741 ms  | 118.675 ms | 9.99 MiB | 21105  |

PlanktonIndividuals.jl now can also run on Apple M-series GPU. Below is a similar benchmark on Apple CPU and GPU.

```julia
PlanktonIndividuals v0.7.5
Julia Version 1.11.7

macOS 26.0.0, Darwin 25.0.0

Toolchain:
- Julia: 1.11.7
- LLVM: 16.0.6

Julia packages: 
- Metal.jl: 1.8.0
- GPUArrays: 11.2.3
- GPUCompiler: 1.6.1
- KernelAbstractions: 0.9.38
- ObjectiveC: 3.4.2
- LLVM: 9.4.2
- LLVMDowngrader_jll: 0.6.0+1

1 device:
- Apple M4 Pro 20 GPU cores (64 GiB Unified Memory)
```

| Arch | N     | min        | median     | mean       | max       | memory     | allocs | samples |
| ---- | ----- | ---------- | ---------- | ---------- | --------- | ---------- | ------ | ------- |
|  CPU |  1024 | 718.333 μs | 771.458 μs | 798.633 μs |  1.006 ms | 467.94 KiB |   3679 |      10 |
|  CPU |  4096 |   2.447 ms |   2.518 ms |   2.565 ms |  2.944 ms | 636.72 KiB |   3679 |      10 |
|  CPU | 16384 |   9.089 ms |   9.216 ms |   9.226 ms |  9.497 ms |   1.28 MiB |   3680 |      10 |
|  CPU | 65536 |  37.136 ms |  37.409 ms |  37.412 ms | 37.770 ms |   3.92 MiB |   3601 |      10 |
|  GPU |  1024 |  70.685 ms |  72.641 ms |  73.513 ms | 79.681 ms |   3.47 MiB | 115344 |      10 |
|  GPU |  4096 |  77.897 ms |  81.464 ms |  81.831 ms | 89.131 ms |   3.48 MiB | 115532 |      10 |
|  GPU | 16384 |  80.004 ms |  81.610 ms |  81.961 ms | 88.787 ms |   3.47 MiB | 115508 |      10 |
|  GPU | 65536 |  81.589 ms |  81.993 ms |  82.922 ms | 91.138 ms |   3.48 MiB | 115611 |      10 |


| Arch |     N |  Ns |        min |     median |       mean |        max |    memory | allocs | samples |
| ---- | ----- | --- | ---------- | ---------- | ---------- | ---------- | --------- | ------ | ------- |
|  CPU |  1024 |  32 |   1.650 ms |   1.727 ms |   1.746 ms |   2.009 ms |  2.20 MiB |   4323 |      10 |
|  CPU |  1024 |  64 |   4.302 ms |   4.335 ms |   4.428 ms |   4.902 ms |  7.53 MiB |   4963 |      10 |
|  CPU |  1024 | 128 |  14.636 ms |  14.975 ms |  15.051 ms |  15.846 ms | 21.81 MiB |   6164 |      10 |
|  CPU |  1024 | 256 |  58.705 ms |  59.051 ms |  59.212 ms |  61.558 ms | 81.09 MiB |   8725 |      10 |
|  CPU |  4096 |  32 |   3.334 ms |   3.425 ms |   3.462 ms |   3.846 ms |  2.36 MiB |   4323 |      10 |
|  CPU |  4096 |  64 |   6.199 ms |   6.277 ms |   6.349 ms |   6.809 ms |  7.70 MiB |   4963 |      10 |
|  CPU |  4096 | 128 |  16.575 ms |  16.851 ms |  16.908 ms |  17.764 ms | 21.98 MiB |   6164 |      10 |
|  CPU |  4096 | 256 |  60.498 ms |  61.014 ms |  61.461 ms |  65.441 ms | 81.26 MiB |   8725 |      10 |
|  CPU | 16384 |  32 |   9.935 ms |  10.054 ms |  10.088 ms |  10.442 ms |  3.02 MiB |   4245 |      10 |
|  CPU | 16384 |  64 |  15.174 ms |  15.359 ms |  15.371 ms |  15.721 ms |  8.36 MiB |   4885 |      10 |
|  CPU | 16384 | 128 |  26.593 ms |  26.860 ms |  26.911 ms |  27.408 ms | 22.64 MiB |   6165 |      10 |
|  CPU | 16384 | 256 |  71.047 ms |  71.621 ms |  71.659 ms |  72.294 ms | 81.92 MiB |   8726 |      10 |
|  CPU | 65536 |  32 |  38.325 ms |  38.712 ms |  38.760 ms |  39.282 ms |  5.66 MiB |   4245 |      10 |
|  CPU | 65536 |  64 |  49.727 ms |  50.438 ms |  50.447 ms |  50.888 ms | 10.99 MiB |   4885 |      10 |
|  CPU | 65536 | 128 |  64.143 ms |  64.990 ms |  64.915 ms |  65.964 ms | 25.28 MiB |   6165 |      10 |
|  CPU | 65536 | 256 | 112.613 ms | 113.154 ms | 113.323 ms | 114.930 ms | 84.55 MiB |   8726 |      10 |
|  GPU |  1024 |  32 |  52.458 ms |  53.845 ms |  54.525 ms |  61.694 ms |  3.58 MiB | 118477 |      10 |
|  GPU |  1024 |  64 |  54.444 ms |  55.915 ms |  57.301 ms |  72.506 ms |  3.93 MiB | 130398 |      10 |
|  GPU |  1024 | 128 |  58.973 ms |  60.532 ms |  62.796 ms |  86.138 ms |  4.85 MiB | 159336 |      10 |
|  GPU |  1024 | 256 |  75.546 ms |  78.564 ms |  83.052 ms | 111.422 ms |  7.59 MiB | 241115 |      10 |
|  GPU |  4096 |  32 |  51.781 ms |  53.934 ms |  55.887 ms |  76.180 ms |  3.58 MiB | 118592 |      10 |
|  GPU |  4096 |  64 |  55.067 ms |  56.103 ms |  56.670 ms |  63.328 ms |  3.93 MiB | 130565 |      10 |
|  GPU |  4096 | 128 |  54.806 ms |  60.063 ms |  60.309 ms |  72.531 ms |  4.85 MiB | 159534 |      10 |
|  GPU |  4096 | 256 |  78.388 ms |  79.156 ms |  82.094 ms | 101.848 ms |  7.59 MiB | 241281 |      10 |
|  GPU | 16384 |  32 |  50.982 ms |  53.887 ms |  54.549 ms |  61.776 ms |  3.58 MiB | 118608 |      10 |
|  GPU | 16384 |  64 |  55.210 ms |  56.055 ms |  58.040 ms |  77.365 ms |  3.93 MiB | 130575 |      10 |
|  GPU | 16384 | 128 |  59.182 ms |  60.448 ms |  62.112 ms |  78.982 ms |  4.85 MiB | 159542 |      10 |
|  GPU | 16384 | 256 |  76.740 ms |  78.906 ms |  81.357 ms | 105.150 ms |  7.59 MiB | 241292 |      10 |
|  GPU | 65536 |  32 |  52.495 ms |  53.323 ms |  54.469 ms |  63.653 ms |  3.58 MiB | 118715 |      10 |
|  GPU | 65536 |  64 |  54.705 ms |  56.284 ms |  57.580 ms |  70.857 ms |  3.93 MiB | 130660 |      10 |
|  GPU | 65536 | 128 |  54.231 ms |  56.561 ms |  57.933 ms |  70.151 ms |  4.85 MiB | 159663 |      10 |
|  GPU | 65536 | 256 |  76.922 ms |  79.406 ms |  82.436 ms | 110.916 ms |  7.60 MiB | 241376 |      10 |



| Arch |     N |  Ns |        min |     median |       mean |        max |     memory |   allocs | samples |
| ---- | ----- | --- | ---------- | ---------- | ---------- | ---------- | ---------- | -------- | ------- |
|  CPU |  1024 |  32 |  13.701 ms |  14.026 ms |  14.027 ms |  14.352 ms |   1.07 MiB |     4311 |      10 |
|  CPU |  1024 |  64 | 104.133 ms | 104.601 ms | 104.670 ms | 105.534 ms |   5.55 MiB |     8023 |      10 |
|  CPU |  1024 | 128 | 828.092 ms | 829.160 ms | 828.977 ms | 829.878 ms |  43.23 MiB |    21591 |       7 |
|  CPU |  1024 | 256 |    6.707 s |    6.707 s |    6.707 s |    6.707 s | 322.40 MiB |    73303 |       1 |
|  CPU |  4096 |  32 |  15.903 ms |  16.124 ms |  16.141 ms |  16.600 ms |   1.24 MiB |     4311 |      10 |
|  CPU |  4096 |  64 | 106.679 ms | 107.173 ms | 107.316 ms | 108.028 ms |   5.72 MiB |     8023 |      10 |
|  CPU |  4096 | 128 | 820.950 ms | 821.295 ms | 821.341 ms | 821.975 ms |  41.12 MiB |    21591 |       7 |
|  CPU |  4096 | 256 |    6.801 s |    6.801 s |    6.801 s |    6.801 s | 322.56 MiB |    73303 |       1 |
|  CPU | 16384 |  32 |  25.392 ms |  25.647 ms |  25.625 ms |  25.994 ms |   1.90 MiB |     4312 |      10 |
|  CPU | 16384 |  64 | 117.064 ms | 117.777 ms | 117.660 ms | 118.296 ms |   6.38 MiB |     8024 |      10 |
|  CPU | 16384 | 128 | 833.843 ms | 835.038 ms | 834.786 ms | 835.743 ms |  41.77 MiB |    21592 |       6 |
|  CPU | 16384 | 256 |    6.712 s |    6.712 s |    6.712 s |    6.712 s | 323.22 MiB |    73304 |       1 |
|  CPU | 65536 |  32 |  65.352 ms |  66.282 ms |  66.218 ms |  66.924 ms |   4.54 MiB |     4312 |      10 |
|  CPU | 65536 |  64 | 165.360 ms | 166.149 ms | 166.157 ms | 167.143 ms |   9.02 MiB |     8024 |      10 |
|  CPU | 65536 | 128 | 888.348 ms | 889.139 ms | 889.124 ms | 890.255 ms |  44.41 MiB |    21592 |       6 |
|  CPU | 65536 | 256 |    6.776 s |    6.776 s |    6.776 s |    6.776 s | 325.86 MiB |    73304 |       1 |
|  GPU |  1024 |  32 |  21.830 ms |  22.871 ms |  23.228 ms |  28.077 ms |   4.56 MiB |   140410 |      10 |
|  GPU |  1024 |  64 |  29.849 ms |  30.973 ms |  32.603 ms |  48.207 ms |  13.65 MiB |   383611 |      10 |
|  GPU |  1024 | 128 | 119.130 ms | 124.133 ms | 126.141 ms | 150.394 ms |  87.53 MiB |  2247365 |      10 |
|  GPU |  1024 | 256 | 932.425 ms | 962.754 ms | 966.198 ms |    1.023 s | 658.84 MiB | 17009358 |       6 |
|  GPU |  4096 |  32 |  21.630 ms |  22.592 ms |  22.975 ms |  26.425 ms |   4.57 MiB |   140587 |      10 |
|  GPU |  4096 |  64 |  29.367 ms |  29.591 ms |  32.406 ms |  54.415 ms |  13.66 MiB |   383807 |      10 |
|  GPU |  4096 | 128 | 118.579 ms | 125.494 ms | 127.510 ms | 154.465 ms |  86.75 MiB |  2247534 |      10 |
|  GPU |  4096 | 256 | 915.231 ms | 950.850 ms | 949.480 ms | 969.102 ms | 658.84 MiB | 17009532 |       6 |
|  GPU | 16384 |  32 |  21.454 ms |  22.460 ms |  22.777 ms |  26.804 ms |   4.57 MiB |   140601 |      10 |
|  GPU | 16384 |  64 |  27.957 ms |  29.587 ms |  31.845 ms |  53.485 ms |  13.65 MiB |   383729 |      10 |
|  GPU | 16384 | 128 | 119.955 ms | 124.230 ms | 127.157 ms | 160.412 ms |  85.97 MiB |  2247538 |      10 |
|  GPU | 16384 | 256 | 904.653 ms | 963.670 ms | 948.708 ms | 966.893 ms | 658.84 MiB | 17009526 |       6 |
|  GPU | 65536 |  32 |  20.600 ms |  23.022 ms |  23.021 ms |  26.879 ms |   4.57 MiB |   140708 |      10 |
|  GPU | 65536 |  64 |  26.643 ms |  29.079 ms |  30.251 ms |  44.419 ms |  13.66 MiB |   383811 |      10 |
|  GPU | 65536 | 128 | 120.069 ms | 121.475 ms | 126.132 ms | 163.317 ms |  87.53 MiB |  2247622 |      10 |
|  GPU | 65536 | 256 | 894.671 ms | 943.257 ms | 938.207 ms | 962.530 ms | 658.84 MiB | 17009612 |       6 |
