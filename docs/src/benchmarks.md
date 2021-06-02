
# [Benchmarks](@id benchmarks)

Here we benchmark the model performance in two `Architecture`s.
The number of individuals used in the benchmark are `(2, 2^5, 2^10, 2^15)`.
And we also use different grid resolutions in 2-Dimensional and 3-Dimensional model setup.

## 0-Dimensional model

This is a benchmark of a simple 0-Dimensional model setup without advection of Eulerian tracers. However, the advection of individuals still take the same amount of time whether the velocity field is provided or not.

```julia
PlanktonIndividuals v0.2.0
Julia Version 1.6.0
Commit f9720dc2eb (2021-03-24 12:55 UTC)
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
|  CPU |     2 |   1.066 ms |   1.231 ms |   2.018 ms | 610.95 KiB |   4385 |
|  CPU |    32 |   1.238 ms |   1.390 ms |   1.977 ms | 610.95 KiB |   4385 |
|  CPU |  1024 |   8.462 ms |   8.663 ms |   9.398 ms | 610.89 KiB |   4381 |
|  CPU | 32768 | 226.577 ms | 231.720 ms | 264.892 ms | 610.89 KiB |   4381 |
|  GPU |     2 |   7.221 ms |   8.179 ms |   9.740 ms |   1.57 MiB |  21480 |
|  GPU |    32 |   7.274 ms |   7.624 ms |   8.229 ms |   1.57 MiB |  21471 |
|  GPU |  1024 |   7.305 ms |   7.647 ms |   8.606 ms |   1.57 MiB |  21471 |
|  GPU | 32768 |   7.871 ms |  12.539 ms |  24.820 ms |   1.57 MiB |  21547 |

## 2-Dimensional model

This is the benchmark of a 2-Dimensional model setup with `(Ns, 1, Ns)` grid cells. Here `Ns = [32, 64, 128]`.

```julia
PlanktonIndividuals v0.2.0
Julia Version 1.6.0
Commit f9720dc2eb (2021-03-24 12:55 UTC)
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
|  CPU |     2 |  32 |   4.996 ms |   5.649 ms |   5.808 ms |  2.83 MiB |   4438 |
|  CPU |     2 |  64 |  14.257 ms |  16.509 ms |  17.855 ms |  8.82 MiB |   4438 |
|  CPU |     2 | 128 |  52.537 ms |  56.669 ms |  72.511 ms | 31.85 MiB |   4438 |
|  CPU |    32 |  32 |   4.884 ms |   5.741 ms |   6.073 ms |  2.83 MiB |   4438 |
|  CPU |    32 |  64 |  14.430 ms |  16.662 ms |  17.849 ms |  8.82 MiB |   4438 |
|  CPU |    32 | 128 |  52.890 ms |  53.086 ms |  53.704 ms | 31.85 MiB |   4438 |
|  CPU |  1024 |  32 |  11.078 ms |  12.091 ms |  13.222 ms |  2.83 MiB |   4434 |
|  CPU |  1024 |  64 |  24.283 ms |  24.721 ms |  26.330 ms |  8.82 MiB |   4438 |
|  CPU |  1024 | 128 |  59.773 ms |  60.021 ms |  60.935 ms | 31.85 MiB |   4438 |
|  CPU | 32768 |  32 | 211.933 ms | 216.529 ms | 241.320 ms |  2.83 MiB |   4434 |
|  CPU | 32768 |  64 | 227.456 ms | 238.156 ms | 285.987 ms |  8.82 MiB |   4438 |
|  CPU | 32768 | 128 | 269.776 ms | 292.217 ms | 329.241 ms | 31.85 MiB |   4438 |
|  GPU |     2 |  32 |   7.423 ms |   9.114 ms |  14.216 ms |  1.60 MiB |  20939 |
|  GPU |     2 |  64 |   7.704 ms |   8.495 ms |  12.629 ms |  1.69 MiB |  20947 |
|  GPU |     2 | 128 |   8.074 ms |   8.476 ms |   9.231 ms |  2.07 MiB |  20947 |
|  GPU |    32 |  32 |   7.121 ms |   7.546 ms |   8.266 ms |  1.60 MiB |  20939 |
|  GPU |    32 |  64 |   7.822 ms |   8.457 ms |  11.789 ms |  1.69 MiB |  20947 |
|  GPU |    32 | 128 |   8.189 ms |   8.508 ms |   9.170 ms |  2.07 MiB |  20947 |
|  GPU |  1024 |  32 |   7.155 ms |   7.678 ms |   8.766 ms |  1.60 MiB |  20939 |
|  GPU |  1024 |  64 |   7.853 ms |   8.174 ms |   9.340 ms |  1.69 MiB |  20947 |
|  GPU |  1024 | 128 |   8.189 ms |   8.525 ms |   9.471 ms |  2.07 MiB |  20947 |
|  GPU | 32768 |  32 |   8.120 ms |  10.474 ms |  18.289 ms |  1.60 MiB |  21079 |
|  GPU | 32768 |  64 |   8.214 ms |  11.399 ms |  21.621 ms |  1.69 MiB |  21083 |
|  GPU | 32768 | 128 |   9.314 ms |  12.102 ms |  20.027 ms |  2.07 MiB |  21083 |

## 3-Dimensional model

This is the benchmark of a 3-Dimensional model setup with `(Ns, Ns, Ns)` grid cells. Here `Ns = [32, 64]`.

```julia
PlanktonIndividuals v0.2.0
Julia Version 1.6.0
Commit f9720dc2eb (2021-03-24 12:55 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, broadwell)
  GPU: Tesla P100-PCIE-12GB
```

| Arch |     N |  Ns |        min |       mean |        max |   memory | allocs |
|------|-------|-----|------------|------------|------------|----------|--------|
|  CPU |     2 |  32 |  39.864 ms |  40.277 ms |  43.291 ms | 1.51 MiB |   4206 |
|  CPU |     2 |  64 | 351.139 ms | 353.212 ms | 361.829 ms | 8.56 MiB |   4207 |
|  CPU |    32 |  32 |  39.855 ms |  40.296 ms |  43.522 ms | 1.51 MiB |   4206 |
|  CPU |    32 |  64 | 376.590 ms | 395.313 ms | 417.473 ms | 8.56 MiB |   4207 |
|  CPU |  1024 |  32 |  49.468 ms |  49.872 ms |  51.700 ms | 1.51 MiB |   4206 |
|  CPU |  1024 |  64 | 366.752 ms | 369.094 ms | 380.962 ms | 8.56 MiB |   4207 |
|  CPU | 32768 |  32 | 260.691 ms | 264.792 ms | 294.048 ms | 1.51 MiB |   4206 |
|  CPU | 32768 |  64 | 583.412 ms | 590.231 ms | 623.753 ms | 8.56 MiB |   4207 |
|  GPU |     2 |  32 |   5.922 ms |   8.108 ms |  15.476 ms | 2.56 MiB |  20019 |
|  GPU |     2 |  64 |  10.199 ms |  13.916 ms |  24.250 ms | 9.61 MiB |  20016 |
|  GPU |    32 |  32 |   5.863 ms |   8.657 ms |  14.639 ms | 2.56 MiB |  20019 |
|  GPU |    32 |  64 |  10.313 ms |  13.910 ms |  24.215 ms | 9.61 MiB |  20016 |
|  GPU |  1024 |  32 |   5.856 ms |   8.798 ms |  13.367 ms | 2.56 MiB |  20019 |
|  GPU |  1024 |  64 |   9.996 ms |  13.269 ms |  20.624 ms | 9.61 MiB |  20016 |
|  GPU | 32768 |  32 |   7.491 ms |  10.591 ms |  18.827 ms | 2.56 MiB |  20155 |
|  GPU | 32768 |  64 |  12.336 ms |  16.536 ms |  20.237 ms | 9.61 MiB |  20156 |
