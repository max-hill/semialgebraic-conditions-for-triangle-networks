# About this repository

This repository contains code implementing methods described in the paper:

[1] _Semialgebraic conditions for distinguishing triangles in phylogenetic
networks_, B. Currie, A. K. Englander, J. A. Esparza-Lozano, E. Gross, M.
Hill, C. Long, D. Olds, K. O'Connor, U. Ranasinghe, and C. Sum. (2026).

# Contents

This repository includes the following code for plots and algebraic
computations:

- [Code for Section (The volume of model intersections)](src/model-size-simulations.jl)

* [Code to replicate Figure 9](src/delta-figure-plot.jl)

* [The computation used the ghost-lineage scenario proof](src/ghost-scenario-calculation.m2)

* [The computation showing the existence of a one-parameter family of solutions](src/one-param-family-computation.m2)

* [Computations accompanying Appendix: The Fourier Transform](src/Fourier-transform-appendix.jl)

* To replicate the figure 2 plot, see [this Desmos link](https://www.desmos.com/calculator/xs5aivcjn7) (also saved as a .json file [here](src/desmos-code-for-figure-2.json))

* Code for replicating the _Rhizoplaca_ example can be found in the following files:

  - [replicating the network](src/rhizoplaca/rhizoplaca-instructions.md)
  - [plotting the networks](src/rhizoplaca/rhizoplaca-plot.jl)
  - [verifying the cutoff](src/rhizoplaca/rhizoplaca-cutoff-computation.jl)


# Project directory structure

The assumed directory structure is as follows:

```
gitrepo
├── data
│   ├── best-rhizoplaca-networks.txt
│   └── ...
├── figures
│   ├── ...
│   └── ...
├── README.md
└── src
    ├── delta-figure-plot.jl
    ├── desmos-code-for-figure-2.json
    ├── Fourier-transform-appendix.jl
    ├── ghost-scenario-calculation.m2
    ├── model-size-simulations.jl
    ├── one-param-family-computation.m2
    └── rhizoplaca
        ├── rhizoplaca-cutoff-computation.jl
        ├── rhizoplaca-instructions.md
        └── rhizoplaca-plot.jl
```

# Software license

Software in this repository is licensed under the GPLv3 License:

> Copyright (c) 2026: Max Hill and coauthors
> 
> This program is free software: you can redistribute it and/or modify it
> under the terms of the GNU General Public License as published by the Free
> Software Foundation, either version 3 of the License, or (at your option)
> any later version.
> 
> This program is distributed in the hope that it will be useful, but WITHOUT
> ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
> FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
> more details.
> 
> You should have received a copy of the GNU General Public License along with
> this program. If not, see <https://www.gnu.org/licenses/>.
