# Biogeochemistry

All the Eulerian tracers are advected and diffused by velocity fields provided by `PlanktonSimulation`.
The advection scheme used in the model is [Third Order Direct Space-Time with Flux Limiting](https://mitgcm.readthedocs.io/en/latest/algorithm/adv-schemes.html#third-order-direct-space-time-with-flux-limiting).

```math
\frac{\partial X}{\partial t} = - \nabla \cdot (\bm{u}X) + \nabla \cdot (\bm{K}\nabla X) + S_X
```

where ``\bm{u}=(u,v,w)``, velocity in physical model, ``\bm{K}`` is the mixing coefficient in physical model, and ``S_X`` is the source and sink term of tracer ``X``.

The source and sinks of each tracer,``S_X``, are different andincluding biological transformations, chemical reactions,and external sources and sinks. 

## Carbon Cycle

```math
S_{DIC} = -\sum_j PS_j\cdot n_j + k_{DOC}\cdot DOC + F_C\\
S_{DOC} = k_{POC} \cdot POC + f_{C,m} \cdot \sum_j ((Bm_j+Cq_j)\cdot n_{j,m}) + f_{C,g} \cdot \sum_j ((Bm_j+Cq_j)\cdot n_{j,g}) - k_{DOC} \cdot DOC\\
S_{POC} = (1-f_{C,m}) \cdot \sum_j ((Bm_j+Cq_j)\cdot n_{j,m}) + (1-f_{C,g}) \cdot \sum_j ((Bm_j+Cq_j)\cdot n_{j,g}) - k_{POC} \cdot POC\\
```

where ``n_j`` is the cell number of species ``j``, ``n_{j,m}`` is the dead cell number of species ``j``, ``n_{j,g}`` is the grazed cell number of sepcies ``j``.

## Nitrogen Cycle

```math
S_{HN4} = -\sum_j VNH4_j\cdot n_j + k_{DON}\cdot DON - k_{nit}\cdot NH4\\
S_{NO3} = -\sum_j VNO3_j\cdot n_j + k_{nit}\cdot NH4\\
S_{DON} = k_{PON} \cdot PON + f_{N,m} \cdot \sum_j ((Bm_j*R_{NC}+Nq_j)\cdot n_{j,m}) + f_{N,g} \cdot \sum_j ((Bm_j*R_{NC}+Nq_j)\cdot n_{j,g}) - k_{DON} \cdot DON\\
S_{PON} = (1-f_{N,m}) \cdot \sum_j ((Bm_j*R_{NC}+Nq_j)\cdot n_{j,m}) + (1-f_{N,g}) \cdot \sum_j ((Bm_j*R_{NC}+Nq_j)\cdot n_{j,g}) - k_{PON} \cdot PON\\
```

## Phosphorus Cycle

```math
S_{PO4} = -\sum_j VPO4_j\cdot n_j + k_{DOP}\cdot DOP\\
S_{DOP} = k_{POP} \cdot POP + f_{P,m} \cdot \sum_j ((Bm_j*R_{PC}+Pq_j)\cdot n_{j,m}) + f_{P,g} \cdot \sum_j ((Bm_j*R_{PC}+Nq_j)\cdot n_{j,g}) - k_{DOP} \cdot DOP\\
S_{POP} = (1-f_{P,m}) \cdot \sum_j ((Bm_j*R_{PC}+Pq_j)\cdot n_{j,m}) + (1-f_{P,g}) \cdot \sum_j ((Bm_j*R_{PC}+Pq_j)\cdot n_{j,g}) - k_{POP} \cdot POP\\
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