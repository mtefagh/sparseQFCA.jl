#-------------------------------------------------------------------------------------------

#=
    Purpose:    Finding blocked reactions in metabolic network
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology - Iran
    Date:       July 2022
=#

#-------------------------------------------------------------------------------------------

module SwiftCC
export swiftCC

using Distributed

@everywhere using GLPK, JuMP, COBREXA, LinearAlgebra, SparseArrays

include("../Data Processing/pre_processing.jl")

using .pre_processing

"""
    swiftCC(ModelObject)

A function that finds blocked reactions for a metabolic network.

# INPUTS

- `ModelObject`:        a newly object of MyModel.

# OPTIONAL INPUTS

-

# OUTPUTS

- `blocked_names`:      Blocked reaction Names.

# EXAMPLES

- Full input/output example
```julia
julia> blocked_reactions = swiftCC(ModelObject)
```

See also: `MyModel`, myModel_Constructor(), 'getTolerance()', `reversibility()`, `homogenization()`

"""

@everywhere function swiftCC(ModelObject::MyModel)

    # Exporting data from ModelObject:

    S = ModelObject.S
    Metabolites = ModelObject.Metabolites
    Reactions = ModelObject.Reactions
    m = ModelObject.m
    n = ModelObject.n
    lb = ModelObject.lb
    ub = ModelObject.ub

    # assigning a small value to Tolerance representing the level of error tolerance:

    Tolerance = getTolerance()

    # Determining the reversibility of a reaction:

    irreversible_reactions_id, reversible_reactions_id = reversibility(lb)

    # Determining the number of irreversible and reversible reactions:

    n_irr = length(irreversible_reactions_id)
    n_rev = length(reversible_reactions_id)

    # Calculating the dimensions of the S matrix:

    row_num, col_num = size(S)

    # Finding irreversible blocked reactions in 1 LP problem:

    lb_u = zeros(n_irr)
    ub_u = ones(n_irr)
    model = Model(GLPK.Optimizer)
    @variable(model, lb[i] <= V[i = 1:n] <= ub[i])
    @variable(model, lb_u[i] <= u[i = 1:n_irr] <= ub_u[i])
    @constraint(model, S * V .== 0)
    objective_function = ones(n_irr)'*u
    @objective(model, Max, objective_function)
    @constraint(model, [i in 1:n_irr], u[i] <= V[irreversible_reactions_id[i]])
    optimize!(model)

    irr_blocked_reactions = []
    for i in range(1,n_irr; step=1)
        if isapprox(value(u[i]), 0.0, atol = Tolerance)
            append!(irr_blocked_reactions, irreversible_reactions_id[i])
        end
    end

    # Determining the number of irreversible blocked reactions:

    irr_blocked_num = length(irr_blocked_reactions)

    # Finding reversible blocked reactions:

    # S_Transpose:

    S_transpose = S'

    # Calculating the dimensions of the S_transpose matrix:

    row_num_trans, col_num_trans = size(S_transpose)

    # Removing irrevesible blocked from Stoichiometric Matrix:

    S_transpose_noIrrBlocked = S_transpose[setdiff(1:end, irr_blocked_reactions), :]

    # Constructing the I_reversible Matrix:

    I_reversible = zeros(n, n_rev)

    rev_id = 1
    for col in eachcol(I_reversible)
        col[reversible_reactions_id[rev_id]] = 1.0
        rev_id = rev_id + 1
    end

    # Removing irrevesible blocked from I_reversible Matrix:

    I_reversible = I_reversible[setdiff(1:end, irr_blocked_reactions), :]

    # Calculating the dimensions of the S_transpose_noIrrBlocked and I_reversible matrix :

    S_trans_row, S_trans_col = size(S_transpose_noIrrBlocked)
    I_row, I_col = size(I_reversible)

    # using Gaussian elimination to find reversible blocked reactions:

    X = S_transpose_noIrrBlocked \ I_reversible

    Sol = (S_transpose_noIrrBlocked*X) - I_reversible

    # Calculating the dimensions of the S_transpose_noIrrBlocked matrix:

    row_sol, col_sol = size(Sol)

    # Finding specific columns in Sol:

    c = 0

    rev_blocked_reactions_col = []
    for col in eachcol(Sol)
        c = c + 1
        if isapprox(norm(col), 0, atol = Tolerance)
            append!(rev_blocked_reactions_col, c)
        end
    end

    rev_blocked_reactions = []
    for i in rev_blocked_reactions_col
        append!(rev_blocked_reactions, reversible_reactions_id[i])
    end

    # Merging reversbile blocked reactions and irreversible blocked reactions:

    blocked_index = []
    blocked_index = union(rev_blocked_reactions, irr_blocked_reactions)

    # Returning a list consist of the Ids of the blocked reactions:

    blocked_index = sort(blocked_index)
    return blocked_index
end

end
