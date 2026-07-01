# About this repository

This repository contains code implementing methods described in the paper:

[1] _Semialgebraic conditions for distinguishing triangles in phylogenetic
networks_, B. Currie, A. K. Englander, J. A. Esparza-Lozano, E. Gross, M.
Hill, C. Long, D. Olds, K. O'Connor, U. Ranasinghe, and C. Sum. (2026).

> Software in this repository is licensed under GPLv3: Copyright (c) Max Hill
> and coauthors.

# Contents

This repository includes the following code for plots and algebraic
computations:

- [Code for Section (The volume of model intersections)](src/model-size-simulations.jl)

* [Code to replicate Figure 9](src/delta-figure-plot.jl)

* [The computation used in the ghost-lineage scenario proof](src/ghost-scenario-calculation.m2)

* [The computation showing the existence of a one-parameter family of solutions](src/one-param-family-computation.m2)

* [Computations accompanying Appendix: The Fourier Transform](src/Fourier-transform-appendix.jl)

* To replicate the plot in Figure 2, see [this Desmos link](https://www.desmos.com/calculator/xs5aivcjn7) (also saved as a .json file [here](src/desmos-code-for-figure-2.json))

* Code for replicating the _Rhizoplaca_ example can be found in the following files:

  - [replicating the network](src/rhizoplaca/rhizoplaca-instructions.md)
  - [plotting the networks](src/rhizoplaca/rhizoplaca-plot.jl)
  - [verifying the cutoff](src/rhizoplaca/rhizoplaca-cutoff-computation.jl)


# Project directory structure

The assumed directory structure is as follows:

```
gitrepo
├── data
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


