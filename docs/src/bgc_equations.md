# Biogeochemistry

All Eulerian tracers are advected, diffused, and affected by sources and sinks:

```math
\frac{\partial X}{\partial t} = - \nabla \cdot (\boldsymbol{u}X) + \nabla \cdot (\boldsymbol{K}\nabla X) + S_X
```

where ``\boldsymbol{u}=(u,v,w)`` is the velocity field provided by a physical model (see [Model Simulation](@ref)), ``\boldsymbol{K}`` is the mixing tensor also from the physical model, and ``S_X`` is the source and sink term for tracer ``X``. 

The source and sinks term, ``S_X``, can be different for each tracer and include biological transformations, chemical reactions, and external sources and sinks as detailed below.

The advection scheme used is [Third Order Direct Space-Time with Flux Limiting](https://mitgcm.readthedocs.io/en/latest/algorithm/adv-schemes.html#third-order-direct-space-time-with-flux-limiting).

## Carbon Cycle

```math
S_{DIC} = -\sum_j PS_j\cdot n_j + k_{DOC}\cdot DOC + F_C
```

```math
\begin{align}
S_{DOC} & = k_{POC} \cdot POC + f_{C,m} \cdot \sum_j ((Bm_j+Cq_j)\cdot n_{j,m}) \nonumber\\
        & \quad 
        + f_{C,g} \cdot \sum_j ((Bm_j+Cq_j)\cdot n_{j,g}) - k_{DOC} \cdot DOC \nonumber
\end{align}
```

```math
\begin{align}
S_{POC} & = (1-f_{C,m}) \cdot \sum_j ((Bm_j+Cq_j)\cdot n_{j,m}) \nonumber\\
        & \quad
        + (1-f_{C,g}) \cdot \sum_j ((Bm_j+Cq_j)\cdot n_{j,g}) - k_{POC} \cdot POC \nonumber
\end{align}
```

where ``n_j`` is the cell number of species ``j``, ``n_{j,m}`` is the dead cell number of species ``j``, ``n_{j,g}`` is the grazed cell number of species ``j``.

## Nitrogen Cycle

```math
\begin{align}
S_{HN4} &= -\sum_j VNH4_j\cdot n_j + k_{DON}\cdot DON - k_{nit}\cdot NH4 \nonumber\\
S_{NO3} &= -\sum_j VNO3_j\cdot n_j + k_{nit}\cdot NH4 \nonumber
\end{align}
```

```math
\begin{align}
S_{DON} & = k_{PON} \cdot PON + f_{N,m} \cdot \sum_j ((Bm_j*R_{NC}+Nq_j)\cdot n_{j,m}) \nonumber\\
        & \quad
        + f_{N,g} \cdot \sum_j ((Bm_j*R_{NC}+Nq_j)\cdot n_{j,g}) - k_{DON} \cdot DON \nonumber
\end{align}
```

```math
\begin{align}
S_{PON} & = (1-f_{N,m}) \cdot \sum_j ((Bm_j*R_{NC}+Nq_j)\cdot n_{j,m}) \nonumber\\
        & \quad
        + (1-f_{N,g}) \cdot \sum_j ((Bm_j*R_{NC}+Nq_j)\cdot n_{j,g}) - k_{PON} \cdot PON \nonumber
\end{align}
```

## Phosphorus Cycle

```math
S_{PO4} = -\sum_j VPO4_j\cdot n_j + k_{DOP}\cdot DOP
```

```math
\begin{align}
S_{DOP} & = k_{POP} \cdot POP + f_{P,m} \cdot \sum_j ((Bm_j*R_{PC}+Pq_j)\cdot n_{j,m}) \nonumber \\
        & \quad
        + f_{P,g} \cdot \sum_j ((Bm_j*R_{PC}+Nq_j)\cdot n_{j,g}) - k_{DOP} \cdot DOP \nonumber
\end{align}
```

```math
\begin{align}
S_{POP} & = (1-f_{P,m}) \cdot \sum_j ((Bm_j*R_{PC}+Pq_j)\cdot n_{j,m}) \nonumber\\
        & \quad
        + (1-f_{P,g}) \cdot \sum_j ((Bm_j*R_{PC}+Pq_j)\cdot n_{j,g}) - k_{POP} \cdot POP \nonumber
\end{align}
```

## Parameters

| Symbol            | Param     | Default | Unit              | Description                       |
|-------------------|-----------|---------|-------------------|-----------------------------------|
| ``k_{DOC}``       | kDOC      | 3.8e-7  | ``s^{-1}``        | Remineralization rate of DOC      |
| ``k_{DON}``       | kDON      | 3.8e-7  | ``s^{-1}``        | Remineralization rate of DON      |
| ``k_{DOP}``       | kDOP      | 3.8e-7  | ``s^{-1}``        | Remineralization rate of DOP      |
| ``k_{POC}``       | kPOC      | 3.8e-7  | ``s^{-1}``        | Remineralization rate of POC      |
| ``k_{PON}``       | kPON      | 3.8e-7  | ``s^{-1}``        | Remineralization rate of PON      |
| ``k_{POP}``       | kPOP      | 3.8e-7  | ``s^{-1}``        | Remineralization rate of POP      |
| ``f_{C,m}``       | mortFracC | 0.5     |                   | Fraction of dead C goes to DOM    |
| ``f_{N,m}``       | mortFracN | 0.5     |                   | Fraction of dead N goes to DOM    |
| ``f_{P,m}``       | mortFracP | 0.5     |                   | Fraction of dead P goes to DOM    |
| ``f_{C,g}``       | grazFracC | 0.5     |                   | Fraction of grazed C goes to DOM  |
| ``f_{N,g}``       | grazFracN | 0.5     |                   | Fraction of grazed N goes to DOM  |
| ``f_{P,g}``       | grazFracP | 0.5     |                   | Fraction of grazed P goes to DOM  |
