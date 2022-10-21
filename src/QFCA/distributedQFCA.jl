#-------------------------------------------------------------------------------------------

#=
    Purpose:    QFCA computes the table of flux coupling relations for a metabolic network
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology - Iran
    Date:       October 2022
=#

#-------------------------------------------------------------------------------------------

module DistributedQFCA
export distributedQFCA

using Distributed

@everywhere using GLPK, JuMP, COBREXA, LinearAlgebra, SparseArrays, Distributed, SharedArrays

include("../Data Processing/pre_processing.jl")

using .pre_processing

include("../Consistency Checking/SwiftCC.jl")

using .SwiftCC

"""
    distributedQFCA(myModel)

Function computes the table of flux coupling relations for a metabolic network.

# INPUTS

- `myModel`:        A model that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

-

# OUTPUTS

- `fctable`:       The resulting flux coupling matrix.
                   The meaning of the entry (i, j) is:
                      * 0 - uncoupled reactions
                      * 1 - fully coupled reactions
                      * 2 - partially coupled reactions
                      * 3 - reaction i is directionally coupled to reaction j
                      * 4 - reaction j is directionally coupled to reaction i

# EXAMPLES

- Full input/output example
```julia
julia> fctable = distributedQFCA(myModel)
```

See also: `dataOfModel()`, `reversibility()`, `homogenization()`, `MyModel`, `myModel_Constructor()`

"""


function QFCA(myModel)

    # Exporting data from StandardModel:

    S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel)

    # Determining the reversibility of a reaction:

    irreversible_reactions_id, reversible_reactions_id = reversibility(lb)

    # Determining the number of irreversible and reversible reactions:

    n_irr = length(irreversible_reactions_id)
    n_rev = length(reversible_reactions_id)

    # Homogenizing the upper_bound and lower_bound of reactions:

    lb, ub = homogenization(lb, ub)

    # Creating a newly object of MyModel:

    ModelObject = MyModel(S, Metabolites, Reactions, Genes, m, n, lb, ub)

    # Removing blocked reactions from metabolic network:

    blocked_index  = swiftCC(ModelObject)
    n_blocked = length(blocked_index)
    Reactions_noBlocked = Reactions[setdiff(range(1, n), blocked_index)]
    lb_noBlocked = lb[setdiff(1:end, blocked_index)]
    ub_noBlocked = ub[setdiff(1:end, blocked_index)]
    S_noBlocked  = S[:, setdiff(1:end, blocked_index)]
    row_noBlocked, col_noBlocked = size(S_noBlocked)

    # Determining Coupling by using swiftCC function:

    @everywhere D = Dict()

    # Defining a SharedMatrix:

    @time begin

    DC_Matrix = SharedArray{Int,2}((col_noBlocked, col_noBlocked), init = false)

    # Finding Coupling between reactions by using swiftCC functions in Parallel:

    @sync @distributed for i in range(1, col_noBlocked)
        lb_temp = lb_noBlocked[i]
        ub_temp = ub_noBlocked[i]
        lb_noBlocked[i] = 0.0
        ub_noBlocked[i] = 0.0
        myModel_Constructor(ModelObject ,S_noBlocked, Metabolites, Reactions_noBlocked, Genes, row_noBlocked, col_noBlocked, lb_noBlocked, ub_noBlocked)
        blocked = swiftCC(ModelObject)
        D[i] = blocked
        D_values = collect(values(D[i]))
        DC_Matrix[i,D_values] .= 1.0
        lb_noBlocked[i] = lb_temp
        ub_noBlocked[i] = ub_temp
    end

    end

    # Defining a matrix to save coupling relations:

    fctable = zeros(col_noBlocked, col_noBlocked)

    # Directional Coupling:

    # i-->j

    for i in range(1, col_noBlocked)
        for j in range(1, col_noBlocked)
            if DC_Matrix[i, j] == 1
                fctable[j,i] = 3.0
            end
        end
    end

    # Self-Coupling:

    fctable[diagind(fctable)] .= 1.0

    # Distingushing between (DC, FC or PC) kinds of coupling:

     for i in range(1, col_noBlocked)
        for j in range(i+1, col_noBlocked)

            # i--->j

            if fctable[i,j] == 3.0 && fctable[j,i] == 0.0
                fctable[j,i] = 4.0

            # j--->i

           elseif fctable[i,j] == 0.0 && fctable[j,i] == 3.0
                fctable[i,j] = 4.0

            # i--->j && j--->i   ----> Partially Coupling

            elseif fctable[i,j] == 3.0 && fctable[j,i] == 3.0
                fctable[i,j] = fctable[j,i] = 2.0
            else
                continue
            end
        end
    end

    # Determining Partial Couples:

    PC = Dict()

    s = 1
    Partially_Couples = []

    for i in range(1, col_noBlocked)
        for j in range(1, col_noBlocked)
            if fctable[i,j] == fctable[j,i] == 2.0
                Partially_Couples = i,j
                PC[s] = Partially_Couples
                s = s + 1
            end
        end
    end

    # Solving a LU for each pair to determine Fully Coupling:

    for key in sort(collect(keys(PC)))

        # Transposing S:

        S_noBlocked_transpose = S_noBlocked'

        # CO

        CO = zeros(col_noBlocked)

        # Setting CO vector:

        CO[PC[key][1]] = -1

        # Removing C:

        S_noBlocked_transpose = S_noBlocked_transpose[1:end .!= PC[key][2], :]

        CO = CO[setdiff(1:end, PC[key][2])]

        # Using Gaussian elimination:

        X = S_noBlocked_transpose \ CO

        Sol = (S_noBlocked_transpose*X) - (CO)

        # Determining Fully Coupling:

        if isapprox(norm(Sol), 0.0, atol = 1e-8)
            fctable[PC[key][1],PC[key][2]] = 1.0
        end
    end

    return fctable

end
