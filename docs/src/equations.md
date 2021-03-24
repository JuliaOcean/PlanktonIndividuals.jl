# parameterization of phytoplankton physiology

![skematic](PI_Quota.jpeg)

## Photosynthesis

Basically, we follow Geider et al. 1998 for the parameterization of photosynthesis, but without nutrient limitation.
**Unit: mmolC/cell/s**  

```math
PP_{max}= PP_{max}^a \cdot Sz^{PP^b}
```

``PP_{max}`` is scaled by a power-law relationship of cell size (``Sz``).

```math
PP=PP_{max}\cdot (1-e^{\frac{-\alpha \cdot I\cdot Chl}{PP_{max}\cdot C}})
PS=PP\cdot Bm
```

## Nutrient Uptake

Include intracellular nutrient limitation:

```math
V_i^{max}= V_i^a \cdot Sz^{V_i^b}
```

```math
RegQ_i=\bigg[\frac{R_{iC}^{max}-Q_i}{R_{iC}^{max}-R_{iC}^{min}}\bigg]_0^1\\
V_i=V_i^{max}\cdot regQ_i\cdot\frac{[i]}{[i]+K_i^{sat}}\cdot Bm
```

**Unit: mmolN/cell/s**, ``i`` denotes ``NH_4``, ``NO_3``, ``PO_4``. ``Q_i`` is the quota of ``i``.
``R_{iC}^{max}`` is the maximum ``i``quota, ``R_{iC}^{min}`` is the minimum ``i``quota.
``K_i^{sat}`` is the half-saturation concentration. ``Bm`` is the cellular functional biomass.

## Biosynthesis, Respiration, and Excretion

### Update reserves

Model will first update C, N, and P quotas used to calculate biosynthesis rate and excretion rate.

```math
Q_C^R = Q_C^R+PP\\
Q_N^R=Q_N^R+V_{NO_3}+V_{NH_4}\\
Q_P^R=Q_P^R+V_{PO_4}\\
```

### Biosynthesis

Potential biosynthesis rates based on C, N, P quotas are calculated below.

```math
k_{mtb}= k_{mtb}^a \cdot Sz^{k_{mtb}^b}
```

```math
BioSynC = k_{mtb}\cdot Q_C^R\\
BioSynN = k_{mtb}\cdot Q_N^R/R_{NC}\\
BioSynP = k_{mtb}\cdot Q_P^R/R_{PC}
```

Then the minimum of these rates are selected as the actual biosynthesis rate.
The difference between carbon quota based biosynthesis rate and actual biosynthesis rate is excretion rate.

```math
BioSyn=min(BioSynC,BioSynN,BioSynP)\\
ExcretC=BioSynC - BioSyn\\
```

### Respiration

```math
Respir = k_{respir}^a \cdot Sz^{k_{respir}^b} \cdot Bm
```

### Update reserves and biomass

```math
Q_C^B = Q_C^B + BioSyn\\
Q_C^R = Q_C^R - BioSyn - Respir\\
Q_N^R = Q_N^R - BioSyn*R_{NC}\\
Q_P^R = Q_N^R - BioSyn*R_{PC}\\
```

### Cell division

We use relative cell size ``RCS`` to indicate cell division.
``RCS`` of the smallest cell is 1.0. ``Q_C^B`` of the smallest cell is ``Cquota``.
Cells start to divide at ``RCS=2.0``. The probability of individual cell division is a sigmoidal function of ``RCS``.

```math
RCS = Q_C^B / Cquota\\
P_{divide} = rand(Bernoulli(0.2*(tanh(a*(RCS-b))+1)))
```

``P_{divide}`` is computed every hour no matter what time step the model is.
