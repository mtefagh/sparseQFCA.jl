#-------------------------------------------------------------------------------------------

#=
    Purpose:    Finding blocked reactions in metabolic network
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology - Iran
    Date:       July 2022
=#

#-------------------------------------------------------------------------------------------

module SwiftConsistencyChecking
export swiftCC

using GLPK, JuMP, COBREXA, LinearAlgebra, SparseArrays

include("../pre_processing.jl")

using .pre_processing

"""
    swiftCC(myModel)

Function that finds blocked reactions in metabolic network.

# INPUTS

- `myModel`:        A model that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

-

# OUTPUTS

- `blocked_names`:  Blocked reaction Names.

# EXAMPLES

- Full input/output example
```julia
julia> blocked_reactions = swiftCC(myModel)
```

See also: `dataOfModel()`, 'getTolerance()', `reversibility()`, `homogenization()`

"""


function swiftCC(myModel)

    # Exporting data from model:

    S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel)

    # assigning a small value to Tolerance representing the level of error tolerance:

    Tolerance = getTolerance()

    # Determining the reversibility of a reaction:

    irreversible_reactions_id, reversible_reactions_id = reversibility(lb)

    # Determining the number of irreversible and reversible reactions:

    n_irr = length(irreversible_reactions_id)
    n_rev = length(reversible_reactions_id)

    # Homogenizing the upper_bound and lower_bound of reactions:

    lb, ub = homogenization(lb, ub)

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
    for i in range(1,n_irr)
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

    unit_vector(i,n) = [zeros(i-1); 1 ; zeros(n-i)]
    reversible_reactions_id = sort(reversible_reactions_id)
    I_reversible = unit_vector(reversible_reactions_id[1], n)
    for i in range(2, n_rev)
        a = unit_vector(reversible_reactions_id[i], n)
        I_reversible = hcat(I_reversible, a)
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

    blocked_names = []
    for i in blocked_index
        r_name = reactions(myModel)[i]
        push!(blocked_names, r_name)
    end

    # Returning a list consist of the names of the blocked reactions:

    blocked_names = sort(blocked_names)
    return blocked_names
end

end
