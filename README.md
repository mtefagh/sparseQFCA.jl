# sparseQFCA
*sparseQFCA* is a Julia package for the sparse [Quantitative Flux Coupling Analysis](https://mtefagh.github.io/qfca/).
The example data files are extracted from the [core *E. coli* model](http://systemsbiology.ucsd.edu/Downloads/EcoliCore).

## Usage
### `certificates, blocked, fctable = QFCA(S, rev)`

### Inputs
* `S`: the associated sparse **stoichiometric matrix**
* `rev`: the boolean vector with trues corresponding to the **reversible reactions**

### Outputs
* `certificates`: the fictitious metabolites for the **sparse positive certificates**
* `blocked`: the boolean vector with trues corresponding to the **blocked reactions**
* `fctable`: the resulting **flux coupling matrix**

## License
sparseQFCA is distributed under the [GNU General Public License v3.0](http://www.gnu.org/copyleft/gpl.html).
