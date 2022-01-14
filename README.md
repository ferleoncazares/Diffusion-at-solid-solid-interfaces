# Diffusion at solid-solid interfaces
Finite volume methods implemented in Matlab to model 1D solute diffusion at solid-solid interfaces, accounting for energy barriers and interfacial segregation. The derivation of the models has been published in: [F.D. León-Cázares and E.I. Galindo-Nava. Phys. Rev. Mat., 5 (2021) 123802](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.123802). Please cite this publication if you benefit from this repository.

[![View Diffusion-at-solid-solid-interfaces on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/105235-diffusion-at-solid-solid-interfaces)

## Requirements
Coded in Matlab R2021a. There is no need for additional packages.

## License
This repository is published under a GNU GPL v3 license ([![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)). See `LICENSE.txt` for full licensing information.

## Contents and how to use

The **Examples.m** file showcases the usage of all the functions in this repository. The user can open the example of their preference and change the input parameters. If the simulation diverges, a smaller convergence criterion ~D\*dt/dx^2<0.5 (printed by the functions) should be selected. Relevant plots of the concentrations and fluxes are also added for each diffusion case.

### Functions - Diffusion simulations 

These functions solve the 1D diffusion equations for a compendium of model bicrystal interfaces. The geometry can be set to cartesian, cylindrical or spherical coordinates. Boundary conditions *'perm'* of constant concentrations at both surfaces were implemented for all these cases:

- **traps1D_perm.m** - Interfacial energy barrier.
- **traps1D_perm_t.m** - Monolayer interfacial trap.
- **traps1D_perm_t_diff.m** - Diffuse interfacial trap.
- **traps1D_perm_t_inh.m** - Monolayer inhomogeneous interfacial trap.

Additionally, the monolayer interfacial trap case is implemented for other two sets of boundary conditions: *'open system'* with zero flux (left) and constant concentration (right) surfaces, and *'closed system'* with two no flux surfaces.

- **traps1D_open_t.m** - Monolayer interfacial trap.
- **traps1D_closed_t.m** - Monolayer interfacial trap.

### Functions - Plotting energy landscapes 

- **energy_landscape_plot.m** - Plotting an energy landscape of a monolayer or no trap.
- **energy_landscape_plot_diff.m** - Plotting the energy landscape of a diffuse trap.
