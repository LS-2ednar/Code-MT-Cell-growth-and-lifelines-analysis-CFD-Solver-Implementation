# Cell growth and lifelines analysis CFD Solver Implementation
---
---
## Content

Examples for the three used examples.

---
---

The following OpenFOAM solvers are provided:
- **plotParticleFoam**
This "solver" is a simplified version of the particleFoam solver which just plots the particles in the system without updating their postions or doing any correction.

- **growthParticleFoam**
Should have been a very complicated solver which should have added functionality to the particleFoam solver so that the growth of new cells in a system should have been modeled, however there where a few diffuclties which ulitmately did not allow the solver to operate as intended. Further, this approach would not be computatonally efficient which is why the concept was canceled. 
 
---
---

The following python scripts are provided:
- **generate_cloudPositions.py**
This script generates a particlePostions file for OpenFOAM based on an initial Grid study case.

- **growth_calculation.py**
Using the long-term steady-state distribuion and the mapping form the other scripts. This scripts determines when and where to add particles to the system based on an stochastic approach

- **particlelifelines.py**
This script tracks the trajectories of particles in the system and maps them to the growth regime mappings to determine the particles distribuiton in the tested system and generates markov chains describing the system, and the long-term steady-state distribution of particles.  

- **mapping.py**
Generates a regime mapping based for hydrodynamic parameters shear stress, normal stress, and kolmogorov length scale. Further, a mapping for the substrate uptake rate is generated. Finally, this script provides two growth regime mappings which average or sum-up the previously described approches.

 ---
