# Phytoplankton Physiology

The various resources and processes are summarized in the schematic below, and further detail is provided in the following sub-sections.

![skematic](PI_Quota.jpeg)

## State Variables

8 state variables are simulated in each phytoplankton individual (see table below)

|Symbol | Unit            | Description                       |
|-------|-----------------|-----------------------------------|
| Sz    |                 | cell size                         |
| Bm    | ``mmol~C/cell`` | functional biomass pool in C      |
| Cq    | ``mmol~C/cell`` | C reserve                         |
| Nq    | ``mmol~N/cell`` | N reserve                         |
| Pq    | ``mmol~P/cell`` | P reserve                         |
| chl   | ``mg~Chl/cell`` | Chla pool                         |
| gen   |                 | generation                        |
| age   |  ``hour``       | age of the individual             |

## Photosynthesis

The parameterization for photosynthesis is shown below.

``PC_{max}`` is scaled by a power-law relationship of cell size (``Sz``).

```math
PC_{max}= PCmax \cdot Sz^{PC_b}
```

```math
PC=PC_{max}\cdot (1-e^{\frac{-\alpha \cdot I\cdot Chl}{PC_{max}\cdot Bm}})
```

```math
PS=PC \cdot Bm
```

``PS`` is cell-specific photosynthesis rate (mmol C/cell/s).

## Nutrient Uptake

The uptakes of various nutrients include intracellular nutrient limitation:

```math
\begin{align}
VNH4_{max} &= VNH4max \cdot Sz^{VN_b} \nonumber \\
VNO3_{max} &= VNO3max \cdot Sz^{VN_b} \nonumber \\
VPO4_{max} &= VPO4max \cdot Sz^{VP_b} \nonumber
\end{align}
```

```math
\begin{align}
Q_N &= (Nq + Bm \cdot R_{NC}) / (Cq + Bm) \nonumber \\
Q_P &= (Pq + Bm \cdot R_{PC}) / (Cq + Bm) \nonumber
\end{align}
```

```math
\begin{align}
regQ_N &= \bigg[\frac{Nqmax-Q_N}{Nqmax - Nqmin}\bigg]_0^1 \nonumber \\
regQ_P &= \bigg[\frac{Pqmax-Q_P}{Pqmax - Pqmin}\bigg]_0^1 \nonumber
\end{align}
```

```math
\begin{align}
VNH4 &= VNH4_{max}\cdot regQ_N\cdot\frac{[NH4]}{[NH4]+K_{NH4}^{sat}}\cdot Bm \nonumber \\
VNO3 &= VNO3_{max}\cdot regQ_N\cdot\frac{[NO3]}{[NO3]+K_{NO3}^{sat}}\cdot Bm \nonumber \\
VPO4 &= VPO4_{max}\cdot regQ_P\cdot\frac{[PO4]}{[PO4]+K_{PO4}^{sat}}\cdot Bm \nonumber
\end{align}
```

``VNH4``, ``VNO3``, and ``VPO4`` are cell-specific uptake rates (mmol N/cell/s or mmol P/cell/s).

## Reserve Update

Model first updates C, N, and P reserves based on photosynthesis rate (``PS``) and nutrient uptake rates (``VNH4``, ``VNO3``, and ``VPO4``). The result is used to calculate the biosynthesis and excretion rates.

```math
\begin{align}
Cq &= Cq+PS \cdot \Delta T \nonumber \\
Nq &= Nq+VNO3+VNH4 \cdot \Delta T \nonumber \\
Pq &= Pq+VPO4 \cdot \Delta T \nonumber
\end{align}
```

## Biosynthesis

Potential biosynthesis rates based on C, N, P quotas are calculated below.

```math
k_{mtb}= kmtb_{max} \cdot Sz^{kmtb_b}
```

```math
\begin{align}
BS_C &= Cq \cdot k_{mtb} \nonumber \\
BS_N &= Nq/R_{NC} \cdot k_{mtb} \nonumber \\
BS_P &= Pq/R_{PC} \cdot k_{mtb} \nonumber
\end{align}
```

```math
BS = min(BS_C, BS_N, BS_P)
```

```math
ExuC = BS_C - BS
```

The minimum of these rates gives the actual biosynthesis rate, ``BS`` (mmol C/cell/s), and the difference between carbon-based biosynthesis rate and ``BS`` gives the excretion rate, ``ExuC`` (mmol C/cell/s).

## Chlorophyll Synthesis

```math
S_{chl} = \rho_{chl} * BS * R_{NC}
```

```math
\begin{equation}
\rho_{chl} =
    \begin{cases}
        chl:N * \frac{PC*Bm}{\alpha I \cdot chl}, & \quad \alpha I > 0,\\
        0 & \quad else.
    \end{cases} \nonumber
\end{equation}
```

## Respiration

```math
Respir = respir_a \cdot Sz^{respir_b} \cdot Bm
```

## Biomass Update

Biosynthesis yields a biomass update and corresponding updates in nutrients. The carbon reserve is further modified by respiration.

```math
\begin{align}
Bm  &= Bm + BS \cdot \Delta T \nonumber \\
Cq  &= Cq - (BS - Respir) \cdot \Delta T \nonumber \\
Nq  &= Nq - BS*R_{NC} \cdot \Delta T \nonumber \\
Pq  &= Pq - BS*R_{PC} \cdot \Delta T \nonumber \\
chl &= chl + S_{chl} \cdot \Delta T \nonumber
\end{align}
```

## Cell division

Cell size, ``Sz``, is used to indicate cell division. For the smallest cell, ``Sz=1.0`` and ``Bm=Cquota``. Cells start to divide at ``Sz=2.0`` and the probability of individual cell division is given by a sigmoidal function of ``Sz``.

```math
Sz = (Bm + Cq) / Cquota
```

```math
P_{divide} = rand(Bernoulli(P_{dvid}*(tanh(stp_{dvid}*(Sz-b))+1)))
```

``P_{divide}`` is computed every hour no matter what time step the model is.

## Parameters

| Symbol            | Param   | Default | Unit              | Description                       |
|-------------------|---------|---------|-------------------|-----------------------------------|
| ``PCmax``         | PCmax   | 4.2e-5  | ``s^{-1}``        | Maximum photosynthesis rate       |
| ``VNH4max``       | VNH4max | 6.9e-6  | ``s^{-1}``        | Maximum ammonium uptake rate      |
| ``VNO3max``       | VNO3max | 6.9e-6  | ``s^{-1}``        | Maximum nitrate uptake rate       |
| ``VPO4max``       | VPO4max | 1.2e-6  | ``s^{-1}``        | Maximum phosphate uptake rate     |
| ``PC_b``          | PC_b    | 0.6     |                   | Shape parameter for PC            |
| ``VN_b``          | VN_b    | 0.6     |                   | Shape parameter for VNH4 and VNO3 |
| ``VP_b``          | VP_b    | 0.6     |                   | Shape parameter for VPO4          |
| ``K^{sat}_{NH4}`` | ksatNH4 | 0.005   | ``mmol~N/m^3``    | Half-saturation constant for NH4  |
| ``K^{sat}_{NO3}`` | ksatNO3 | 0.010   | ``mmol~N/m^3``    | Half-saturation constant for NO3  |
| ``K^{sat}_{PO4}`` | ksatPO4 | 0.003   | ``mmol~P/m^3``    | Half-saturation constant for PO4  |
| ``Nqmax``         | Nqmax   | 0.12    | ``mmol~N/mmol~C`` | Maximum N quota in cell           |
| ``Nqmin``         | Nqmin   | 0.05    | ``mmol~N/mmol~C`` | Minimum N quota in cell           |
| ``Pqmax``         | Pqmax   | 0.01    | ``mmol~P/mmol~C`` | Maximum P quota in cell           |
| ``Pqmin``         | Pqmax   | 0.004   | ``mmol~P/mmol~C`` | Minimum P quota in cell           |
| ``R_{NC}``        | R_NC    | 16/106  | ``mmol~N/mmol~C`` | N:C ratio in function biomass     |
| ``R_{PC}``        | R_PC    | 1/106   | ``mmol~P/mmol~C`` | P:C ratio in function biomass     |
| ``kmtb_{max}``    | k_mtb   | 3.5e-5  | ``s^{-1}``        | Maximum metabolic rate            |
| ``kmtb_b``        |k\_mtb\_b| 0.25    |                   | Shape parameter for k_mtb         |
| ``respir_a``      | respir_a| 1.2e-6  | ``s^{-1}``        | Maximum respiration rate          |
| ``respir_b``      | respir_b| 0.6     |                   | Shape parameter for respir_a      |
| ``chl\text{:}N``  | Chl2N   | 3.0     | ``mg~chl/mmol~N`` | Maximum Chl:N in cell             |
| ``P_{dvid}``      | P_dvid  | 5.0e-5  | ``s^{-1}``        | Probability of division per second|
| ``stp_{dvid}``    | stp_dvid| 6.0     |                   | Steepness of division function    |
