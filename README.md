# Sparse Quantitative Flux Coupling Analysis

| **Tests** | **Coverage** |
|:---:|:---:|
| [![Build Status](https://api.travis-ci.com/mtefagh/sparseQFCA.jl.svg?branch=master)](https://app.travis-ci.com/mtefagh/sparseQFCA.jl) | [![Coverage Status](https://coveralls.io/repos/github/mtefagh/sparseQFCA.jl/badge.svg?branch=master)](https://coveralls.io/github/mtefagh/sparseQFCA.jl?branch=master) [![codecov](https://codecov.io/gh/mtefagh/sparseQFCA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mtefagh/sparseQFCA.jl)

## sparseQFCA architecture
The *sparseQFCA* repository offers a [<img src="https://julialang.org/assets/infra/logo.svg" height="20" />](https://julialang.org/) is a Julia package for applying convex optimization techniques to constraint-based modeling of metabolic networks. It includes functions for consistency checking, reversibility correction, and flux coupling analysis. The repository also introduces [QuantomeRedNet](https://github.com/mtefagh/sparseQFCA.jl/blob/master/example/QuantomeRedNet.ipynb), a module designed to provide a loss-less method for the quantitative reduction of metabolic networks, which is particularly beneficial for metabolic engineering and can be effectively used in strain design algorithms. All optimization algorithms are modeled using the [JuMP.jl](https://github.com/jump-dev/JuMP.jl) package and are implemented with parallel processing, utilizing [Distributed.jl](https://github.com/JuliaLang/Distributed.jl) for efficient execution.

![sparseQFCA architecture](/example/sparseQFCA.png "sparseQFCA architecture")

## How to get started

### Prerequisites and requirements

- **Operating system**: Use Linux (Debian, Ubuntu or centOS), MacOS, or Windows
  10 as your operating system. `sparseQFCA` has been tested on these systems.
- **Julia language**: In order to use `sparseQFCA`, you need to install Julia 1.0
  or higher. Download and follow the installation instructions for Julia
  [here](https://julialang.org/downloads/).
- **Hardware requirements**: `sparseQFCA` runs on any hardware that can run Julia,
  and can easily use resources from multiple computers interconnected on a
  network. For processing large datasets, you are required to ensure that the
  total amount of available RAM on all involved computers is larger than the
  data size.
- **Optimization solvers**: `sparseQFCA` uses
  [`JuMP.jl`](https://github.com/jump-dev/JuMP.jl) to formulate optimization
  problems and is compatible with all [`JuMP` supported
  solvers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers).
  However, to perform analysis at least one of these solvers needs to be
  installed on your machine. For a pure Julia implementation, you may use e.g.
  [`HiGHS.jl`](https://github.com/jump-dev/HiGHS.jl), but other solvers ([`GLPK.jl`](https://github.com/jump-dev/GLPK.jl),
  [`CPLEX.jl`](https://github.com/jump-dev/CPLEX.jl), [`MosekTools.jl`](https://github.com/jump-dev/MosekTools.jl),...) work just as well.

:bulb: If you are new to Julia, it is advisable to [familiarize yourself with
the environment
first](https://docs.julialang.org/en/v1/manual/getting-started/).  Use the
Julia [documentation](https://docs.julialang.org) to solve various
language-related issues, and the [Julia package manager
docs](https://julialang.github.io/Pkg.jl/v1/getting-started/) to solve
installation-related difficulties. Of course, [the Julia
channel](https://discourse.julialang.org/) is another fast and easy way to find
answers to Julia specific questions.

## Quick Start
To get started, first run `import Pkg; Pkg.add("sparseQFCA")` to install the *sparseQFCA* package.

## License
*sparseQFCA* is distributed under the [GNU General Public License v3.0](http://www.gnu.org/copyleft/gpl.html).
