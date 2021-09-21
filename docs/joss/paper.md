---
title: 'PlanktonIndividuals.jl: A GPU supported individual-based phytoplankton life cycle model.'
tags:
  - Julia
  - Individual-based model
  - Lagrangian particles
  - phytoplankton
  - biogeochemistry
  - GPU support
authors:
  - name: Zhen Wu
    orcid: 0000-0001-8474-4274
    affiliation: "1"
  - name: Gael Forget
    orcid: 0000-0002-4234-056X
    affiliation: "1"

affiliations:
 - name: MIT, EAPS
   index: 1

date: 15 September 2021
bibliography: paper.bib

---

# Summary
PlanktonIndividuals.jl is a fast individual-based phytoplankton life cycle model written in Julia that runs on both CPU and GPU. It simulates the life cycle of ocean phytoplankton cells as Lagrangian particles while nutrients are represented as Eulerian tracers and advected over the gridded domain. 

The model is used to simulate the temporal and spatial variations in phytoplankton cell density, stoichiometry, as well as growth and division behaviors induced by diel cycle and physical motions ranging from sub-mesoscale to large scale processes. The phytoplankton physiology is well resolved in the model (FIGURE) with the widely-used Droop model(CITATION) implemented for nutrient uptakes. The photosynthesis formulation by Geider et al(CITATION) is also implemented for carbon fixation. Additionally, exudation and mixotrophy are also included in the model(CITATION).

![Schematic diagram of phytoplankton physiology described in PlanktonIndividuals.jl.\label{fig:phyto}](PI_Quota.jpeg)

PlanktonIndividuals.jl can simulate multiple functional groups of phytoplankton with different growth and division strategies which will illustrate the interactions within and between functional groups. The individuals can be simulated not only in a zero-dimensional domain (like lab experiments) but also in one-, two- or three-dimensional domains where individuals will be advected by velocities provided by various models or observations (e.g. MITgcm, Oceananigans.jl).


# Statement of need

