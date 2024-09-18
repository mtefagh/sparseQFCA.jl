# Sparse Quantitative Flux Coupling Analysis

| **Tests** | **Coverage** |
|:---:|:---:|
| [![Build Status](https://app.travis-ci.com/mtefagh/sparseQFCA.jl.svg&branch=master) | [![Coverage Status](https://coveralls.io/repos/github/mtefagh/sparseQFCA.jl/badge.svg?branch=master)](https://coveralls.io/github/mtefagh/sparseQFCA.jl?branch=master) [![codecov](https://codecov.io/gh/mtefagh/sparseQFCA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mtefagh/sparseQFCA.jl) 

## sparseQFCA architecture
*sparseQFCA* is a registered [<img src="https://julialang.org/assets/infra/logo.svg" height="20" />](https://julialang.org/) is a Julia package providing a sparse Quantitative flux coupling analysis(QFCA). It uses parallel processing and is designed for finding flux coupling table and metabolic network reductions, specifically for QuantomeRedNet. Moreover, a Julia implementation of Swift Consistency Checking is also available as a preprocessing subroutine.

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
  [`GLPK.jl`](https://www.gnu.org/software/glpk/), but other solvers (Tulip,
  Gurobi, ...) work just as well.

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
