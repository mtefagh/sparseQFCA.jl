# Sparse Quantitative Flux Coupling Analysis
[![Build Status](https://travis-ci.com/mtefagh/sparseQFCA.jl.svg?branch=master)](https://app.travis-ci.com/mtefagh/sparseQFCA.jl)
[![Coverage Status](https://coveralls.io/repos/github/mtefagh/sparseQFCA.jl/badge.svg?branch=master)](https://coveralls.io/github/mtefagh/sparseQFCA.jl?branch=master)
[![codecov](https://codecov.io/gh/mtefagh/sparseQFCA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mtefagh/sparseQFCA.jl)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mtefagh/sparseQFCA.jl/master?labpath=sparseQFCA.jl%2Fexample%2FdistributedQFCA.ipynb)

*sparseQFCA* is a registered [<img src="https://julialang.org/assets/infra/logo.svg" height="20" />](https://julialang.org/) package for the sparse [Quantitative Flux Coupling Analysis](https://mtefagh.github.io/qfca/). Moreover, a Julia implementation of [Swift Consistency Checking](https://mtefagh.github.io/swiftcore/) is also available as a preprocessing subroutine. The example data files `S.csv` and `rev.csv` are extracted from the [core *E. coli* model](http://systemsbiology.ucsd.edu/Downloads/EcoliCore).

## Usage
### `certificates, blocked, fctable = QFCA(S, rev)`

### Inputs
* `S`: the associated sparse **stoichiometric matrix**
* `rev`: the boolean vector with trues corresponding to the **reversible reactions**

### Outputs
* `certificates`: the fictitious metabolites for the sparse **positive certificates**
* `blocked`: the boolean vector with trues corresponding to the **blocked reactions**
* `fctable`: the resulting **flux coupling** matrix

## Quick Start
To get started, first run `import Pkg; Pkg.add("sparseQFCA")` to install the *sparseQFCA* package, and then see [this Jupyter notebook](https://nbviewer.org/github/mtefagh/sparseQFCA/blob/master/example/sparseQFCA.ipynb) for a demonstration on how to use it.

## License
*sparseQFCA* is distributed under the [GNU General Public License v3.0](http://www.gnu.org/copyleft/gpl.html).
