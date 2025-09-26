# [Benchmarks](@id benchmarks)

Here we benchmark the model performance in two `Architecture`s.
The number of individuals used in the benchmark are `(2^10, 2^12, 2^14, 2^15)`.
And we also use different grid resolutions in 2-Dimensional and 3-Dimensional model setup.

## 0-Dimensional model

This is a benchmark of a simple 0-Dimensional model setup without advection of Eulerian tracers. However, the advection of individuals still take the same amount of time whether the velocity field is provided or not.

```julia
PlanktonIndividuals v0.7.5
Julia Version 1.11.7
Commit f2b3dbda30a (2025-09-08 12:10 UTC)
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 56 × Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
  WORD_SIZE: 64
  LLVM: libLLVM-16.0.6 (ORCJIT, skylake-avx512)
  GPU: Quadro GV100 (sm_70, 32.000 GiB available)
  CUDA runtime 12.9, artifact installation
  CUDA driver 565.57.1 for 12.7
```

| Arch |     N |       min |    median |      mean |       max |     memory | allocs | samples |
|------|-------|-----------|-----------|-----------|-----------|------------|--------|---------|
|  CPU |  1024 |  2.536 ms |  2.629 ms |  2.696 ms |  3.204 ms | 463.71 KiB |   3673 |      10 |
|  CPU |  4096 |  8.091 ms |  8.201 ms |  8.252 ms |  8.829 ms | 632.18 KiB |   3673 |      10 |
|  CPU | 16384 | 30.433 ms | 30.558 ms | 30.745 ms | 31.809 ms |   1.28 MiB |   3595 |      10 |
|  CPU | 32768 | 59.959 ms | 60.364 ms | 60.354 ms | 60.980 ms |   2.15 MiB |   3595 |      10 |
|  GPU |  1024 | 13.006 ms | 13.194 ms | 13.322 ms | 14.415 ms |   2.68 MiB |  77257 |      10 |
|  GPU |  4096 | 13.152 ms | 13.334 ms | 13.386 ms | 13.991 ms |   2.68 MiB |  77257 |      10 |
|  GPU | 16384 | 13.562 ms | 13.755 ms | 13.800 ms | 14.595 ms |   2.68 MiB |  77260 |      10 |
|  GPU | 32768 | 14.646 ms | 14.879 ms | 14.948 ms | 15.450 ms |   2.68 MiB |  77263 |      10 |

## 2-Dimensional model

This is the benchmark of a 2-Dimensional model setup with `(Ns, 1, Ns)` grid cells. Here `Ns = [32, 64, 128]`.

```julia
PlanktonIndividuals v0.7.5
Julia Version 1.11.7
Commit f2b3dbda30a (2025-09-08 12:10 UTC)
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 56 × Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
  WORD_SIZE: 64
  LLVM: libLLVM-16.0.6 (ORCJIT, skylake-avx512)
  GPU: Quadro GV100 (sm_70, 32.000 GiB available)
  CUDA runtime 12.9, artifact installation
  CUDA driver 565.57.1 for 12.7
```

| Arch |     N |  Ns |        min |     median |       mean |        max |    memory | allocs | samples |
|------|-------|-----|------------|------------|------------|------------|-----------|--------|---------|
|  CPU |  1024 |  32 |   5.659 ms |   5.788 ms |   5.874 ms |   6.463 ms |  1.92 MiB |   4317 |      10 |
|  CPU |  1024 |  64 |  13.497 ms |  13.608 ms |  13.770 ms |  15.069 ms |  5.80 MiB |   4878 |      10 |
|  CPU |  1024 | 128 |  45.442 ms |  54.734 ms |  53.155 ms |  61.680 ms | 20.65 MiB |   6158 |      10 |
|  CPU |  4096 |  32 |  11.307 ms |  11.427 ms |  11.519 ms |  12.184 ms |  2.08 MiB |   4238 |      10 |
|  CPU |  4096 |  64 |  19.120 ms |  19.479 ms |  19.667 ms |  20.766 ms |  5.96 MiB |   4878 |      10 |
|  CPU |  4096 | 128 |  51.716 ms |  56.439 ms |  57.364 ms |  65.370 ms | 20.82 MiB |   6158 |      10 |
|  CPU | 16384 |  32 |  33.569 ms |  33.907 ms |  34.243 ms |  35.965 ms |  2.74 MiB |   4239 |      10 |
|  CPU | 16384 |  64 |  41.597 ms |  42.110 ms |  42.844 ms |  45.371 ms |  6.62 MiB |   4879 |      10 |
|  CPU | 16384 | 128 |  75.032 ms |  87.652 ms |  83.260 ms |  89.017 ms | 21.47 MiB |   6159 |      10 |
|  CPU | 32768 |  32 |  63.176 ms |  63.657 ms |  63.717 ms |  64.465 ms |  3.62 MiB |   4239 |      10 |
|  CPU | 32768 |  64 |  71.786 ms |  72.317 ms |  73.374 ms |  76.535 ms |  7.50 MiB |   4879 |      10 |
|  CPU | 32768 | 128 | 106.093 ms | 116.530 ms | 113.994 ms | 120.500 ms | 22.35 MiB |   6159 |      10 |
|  GPU |  1024 |  32 |  12.915 ms |  13.093 ms |  13.219 ms |  13.851 ms |  2.87 MiB |  83564 |      10 |
|  GPU |  1024 |  64 |  13.688 ms |  14.272 ms |  14.355 ms |  15.455 ms |  3.15 MiB |  93423 |      10 |
|  GPU |  1024 | 128 |  15.540 ms |  16.113 ms |  16.036 ms |  16.361 ms |  3.92 MiB | 117894 |      10 |
|  GPU |  4096 |  32 |  12.888 ms |  13.249 ms |  13.541 ms |  15.147 ms |  2.87 MiB |  83564 |      10 |
|  GPU |  4096 |  64 |  13.768 ms |  13.912 ms |  14.277 ms |  15.145 ms |  3.15 MiB |  93423 |      10 |
|  GPU |  4096 | 128 |  15.740 ms |  16.697 ms |  16.592 ms |  17.422 ms |  3.92 MiB | 117894 |      10 |
|  GPU | 16384 |  32 |  13.514 ms |  13.832 ms |  14.139 ms |  16.544 ms |  2.87 MiB |  83565 |      10 |
|  GPU | 16384 |  64 |  13.956 ms |  14.564 ms |  14.831 ms |  17.842 ms |  3.15 MiB |  93425 |      10 |
|  GPU | 16384 | 128 |  15.721 ms |  15.844 ms |  15.941 ms |  16.880 ms |  3.92 MiB | 117896 |      10 |
|  GPU | 32768 |  32 |  13.689 ms |  13.823 ms |  13.926 ms |  15.010 ms |  2.87 MiB |  83568 |      10 |
|  GPU | 32768 |  64 |  14.460 ms |  15.067 ms |  15.092 ms |  15.700 ms |  3.15 MiB |  93428 |      10 |
|  GPU | 32768 | 128 |  16.284 ms |  17.469 ms |  17.356 ms |  18.070 ms |  3.92 MiB | 117899 |      10 |

## 3-Dimensional model

This is the benchmark of a 3-Dimensional model setup with `(Ns, Ns, Ns)` grid cells. Here `Ns = [32, 64, 128]`.

```julia
PlanktonIndividuals v0.7.5
Julia Version 1.11.7
Commit f2b3dbda30a (2025-09-08 12:10 UTC)
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 56 × Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz
  WORD_SIZE: 64
  LLVM: libLLVM-16.0.6 (ORCJIT, skylake-avx512)
  GPU: Quadro GV100 (sm_70, 32 GiB available)
  CUDA runtime 12.9, artifact installation
  CUDA driver 565.57.1 for 12.7
```

| Arch |     N |  Ns |        min |     median |       mean |        max |    memory |  allocs | samples |
|------|-------|-----|------------|------------|------------|------------|-----------|---------|---------|
|  CPU |  1024 |  32 |  51.587 ms |  52.013 ms |  52.390 ms |  54.706 ms |  1.06 MiB |    4305 |      10 |
|  CPU |  1024 |  64 | 407.569 ms | 414.229 ms | 413.787 ms | 419.221 ms |  5.55 MiB |    8017 |      10 |
|  CPU |  1024 | 128 |    3.293 s |    3.301 s |    3.301 s |    3.310 s | 40.95 MiB |   21585 |       2 |
|  CPU |  4096 |  32 |  56.841 ms |  58.601 ms |  58.606 ms |  60.820 ms |  1.22 MiB |    4305 |      10 |
|  CPU |  4096 |  64 | 422.655 ms | 425.630 ms | 425.657 ms | 428.760 ms |  5.72 MiB |    8017 |      10 |
|  CPU |  4096 | 128 |    3.318 s |    3.337 s |    3.337 s |    3.357 s | 41.11 MiB |   21585 |       2 |
|  CPU | 16384 |  32 |  81.034 ms |  82.007 ms |  82.366 ms |  84.177 ms |  1.88 MiB |    4306 |      10 |
|  CPU | 16384 |  64 | 454.705 ms | 457.783 ms | 457.526 ms | 459.801 ms |  6.37 MiB |    8018 |      10 |
|  CPU | 16384 | 128 |    3.336 s |    3.374 s |    3.374 s |    3.413 s | 41.77 MiB |   21586 |       2 |
|  CPU | 32768 |  32 | 112.225 ms | 113.137 ms | 113.290 ms | 114.879 ms |  2.76 MiB |    4306 |      10 |
|  CPU | 32768 |  64 | 495.213 ms | 497.812 ms | 498.189 ms | 501.221 ms |  7.25 MiB |    8018 |      10 |
|  CPU | 32768 | 128 |    3.472 s |    3.482 s |    3.482 s |    3.493 s | 42.65 MiB |   21586 |       2 |
|  GPU |  1024 |  32 |  12.288 ms |  12.387 ms |  12.577 ms |  13.428 ms |  4.08 MiB |  115165 |      10 |
|  GPU |  1024 |  64 |  20.906 ms |  21.227 ms |  21.697 ms |  24.127 ms | 13.06 MiB |  353213 |      10 |
|  GPU |  1024 | 128 |  90.657 ms | 110.816 ms | 113.920 ms | 168.241 ms | 83.72 MiB | 2211802 |      10 |
|  GPU |  4096 |  32 |  12.225 ms |  12.322 ms |  12.434 ms |  13.159 ms |  4.08 MiB |  115165 |      10 |
|  GPU |  4096 |  64 |  20.767 ms |  21.059 ms |  21.390 ms |  24.283 ms | 13.06 MiB |  353213 |      10 |
|  GPU |  4096 | 128 |  90.616 ms | 110.673 ms | 114.507 ms | 170.346 ms | 83.72 MiB | 2211802 |      10 |
|  GPU | 16384 |  32 |  12.400 ms |  12.531 ms |  12.670 ms |  13.325 ms |  4.08 MiB |  115167 |      10 |
|  GPU | 16384 |  64 |  21.132 ms |  21.983 ms |  22.558 ms |  25.184 ms | 13.06 MiB |  353215 |      10 |
|  GPU | 16384 | 128 |  90.849 ms | 110.157 ms | 114.308 ms | 169.462 ms | 83.72 MiB | 2211804 |      10 |
|  GPU | 32768 |  32 |  13.077 ms |  14.088 ms |  13.939 ms |  14.860 ms |  4.08 MiB |  115170 |      10 |
|  GPU | 32768 |  64 |  21.697 ms |  22.645 ms |  22.986 ms |  25.235 ms | 13.06 MiB |  353217 |      10 |
|  GPU | 32768 | 128 |  92.342 ms | 110.662 ms | 112.529 ms | 169.928 ms | 83.72 MiB | 2211807 |      10 |

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
