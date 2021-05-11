# Phytoplankton Physiology

The various resources and processes are summarized in the schematic below, and further detail is provided in the following sub-sections.

![skematic](PI_Quota.jpeg)

## Photosynthesis

The parameterization for photosynthesis follows Geider et al. 1998, but without nutrient limitation (which is treated separately). 

```math
PP_{max}= PP_{max}^a \cdot Sz^{PP^b}
```

``PP_{max}`` is scaled by a power-law relationship of cell size (``Sz``).

```math
PP=PP_{max}\cdot (1-e^{\frac{-\alpha \cdot I\cdot Chl}{PP_{max}\cdot C}})
PS=PP\cdot Bm
```

**Units: mmolC/cell/s**  

## Nutrient Uptake

Include intracellular nutrient limitation:

```math
V_i^{max}= V_i^a \cdot Sz^{V_i^b}
```

```math
RegQ_i=\bigg[\frac{R_{iC}^{max}-Q_i}{R_{iC}^{max}-R_{iC}^{min}}\bigg]_0^1\\
V_i=V_i^{max}\cdot regQ_i\cdot\frac{[i]}{[i]+K_i^{sat}}\cdot Bm
```

where ``i`` denotes one nutrient (``NH_4``, ``NO_3``, or ``PO_4``), ``Q_i`` is the cell quota for ``i``, ``R_{iC}^{max}`` is the maximum ``i`` quota, ``R_{iC}^{min}`` is the minimum ``i`` quota, ``K_i^{sat}`` is the half-saturation concentration, and ``Bm`` is the cellular functional biomass.

**Units: mmolN/cell/s**

## Reserve Update

Model first updates C, N, and P quotas based on `PP` and the `V_i` terms. The result is used to calculate the biosynthesis and excretion rates.

```math
Q_C^R=Q_C^R+PP\\
Q_N^R=Q_N^R+V_{NO_3}+V_{NH_4}\\
Q_P^R=Q_P^R+V_{PO_4}\\
```

## Biosynthesis

Potential biosynthesis rates based on C, N, P quotas are calculated as follows.

```math
k_{mtb}= k_{mtb}^a \cdot Sz^{k_{mtb}^b}
```

```math
BioSynC = k_{mtb}\cdot Q_C^R\\
BioSynN = k_{mtb}\cdot Q_N^R/R_{NC}\\
BioSynP = k_{mtb}\cdot Q_P^R/R_{PC}
```

The minimum of these rates gives the actual biosynthesis rate, `BioSyn`, and the difference between carbon-based biosynthesis rate and `BioSyn` gives the excretion rate, `ExcretC`.

```math
BioSyn=min(BioSynC,BioSynN,BioSynP)\\
ExcretC=BioSynC - BioSyn\\
```

## Respiration

```math
Respir = k_{respir}^a \cdot Sz^{k_{respir}^b} \cdot Bm
```

## Biomass Update

Biosynthesis yields a biomass update and corresponding updates in nutrients. The carbon reserve is further modified by respiration.

```math
Q_C^B = Q_C^B + BioSyn\\
Q_C^R = Q_C^R - BioSyn - Respir\\
Q_N^R = Q_N^R - BioSyn*R_{NC}\\
Q_P^R = Q_N^R - BioSyn*R_{PC}\\
```

## Cell division

Relative cell size, `RCS`, is used to indicate cell division. For the smallest cell, ``RCS=1.0`` and ``Q_C^B=Cquota``. Cells start to divide at ``RCS=2.0`` and the probability of individual cell division is given by a sigmoidal function of `RCS`.

```math
RCS = Q_C^B / Cquota\\
P_{divide} = rand(Bernoulli(0.2*(tanh(a*(RCS-b))+1)))
```

``P_{divide}`` is computed every hour no matter what time step the model is.

## References

- Geider et al 1998
