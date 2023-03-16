#-------------------------------------------------------------------------------------------
#=
    Purpose:    Finding blocked reactions in metabolic network
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology - Iran
    Date:       July 2022
=#
#-------------------------------------------------------------------------------------------

module SwiftCC
export MyModel, myModel_Constructor, swiftCC

using GLPK, JuMP, COBREXA, LinearAlgebra, SparseArrays, Distributed

include("../Data Processing/pre_processing.jl")

using .pre_processing


"""
    MyModel(S, Metabolites, Reactions, Genes, m, n, lb, ub)

A general type for storing a StandardModel which contains the following fields:

- `S`:              LHS matrix (m x n)
- `Metabolites`:    List of metabolic network metabolites.
- `Reactions`:      List of metabolic network reactions.
- `Genes`:          List of metabolic network reactions.
- `m`:              Number of rows of stoichiometric matrix.
- `n`:              Number of columns of stoichiometric matrix.
- `lb`:             Lower bound vector (n x 1)
- `ub`:             Upper bound vector (n x 1)

"""

mutable struct MyModel
    S              ::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}
    Metabolites    ::Array{String,1}
    Reactions      ::Array{String,1}
    Genes          ::Array{String,1}
    m              ::Int
    n              ::Int
    lb             ::Array{Float64,1}
    ub             ::Array{Float64,1}
end

"""
    myModel_Constructor(ModelObject, S, Metabolites, Reactions, Genes, m, n, lb, ub)

A function that initializes a newly created object of MyModel.

# INPUTS

-'ModelObject'      a newly object of MyModel.
- `S`:              LHS matrix (m x n)
- `Metabolites`:    List of metabolic network metabolites.
- `Reactions`:      List of metabolic network reactions.
- `Genes`:          List of metabolic network reactions.
- `m`:              Number of rows of stoichiometric matrix.
- `n`:              Number of columns of stoichiometric matrix.
- `lb`:             Lower bound vector (n x 1)
- `ub`:             Upper bound vector (n x 1)

"""

function myModel_Constructor(ModelObject::MyModel, S::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}, Metabolites::Array{String,1}, Reactions::Array{String,1},
                                         Genes::Array{String,1}, m::Int, n::Int, lb::Array{Float64,1}, ub::Array{Float64,1})
     ModelObject.S = S
     ModelObject.Metabolites = Metabolites
     ModelObject.Reactions = Reactions
     ModelObject.Genes = Genes
     ModelObject.m = m
     ModelObject.n = n
     ModelObject.lb = lb
     ModelObject.ub = ub
end
"""
    swiftCC(ModelObject)

A function that finds blocked reactions for a metabolic network.

# INPUTS

- `ModelObject`:        A newly object of MyModel.

# OPTIONAL INPUTS

- `Tolerance`:          A small number that represents the level of error tolerance.

# OUTPUTS

- `blocked_index`:      Blocked reaction Ids.
- `dualVar`:            Dual variables of a specific constraint.

# EXAMPLES

- Full input/output example
```julia
julia> blocked_index, dual_var = swiftCC(ModelObject)
```

See also: `MyModel`, myModel_Constructor(), `reversibility()`, `homogenization()`

"""

function swiftCC(ModelObject::MyModel, Tolerance::Float64=1e-6)

    # Exporting data from ModelObject:

    S = ModelObject.S
    Metabolites = ModelObject.Metabolites
    Reactions = ModelObject.Reactions
    m = ModelObject.m
    n = ModelObject.n
    lb = ModelObject.lb
    ub = ModelObject.ub

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
    @constraint(model, c1, S * V .== 0)
    objective_function = ones(n_irr)'*u
    @objective(model, Max, objective_function)
    @constraint(model, [i in 1:n_irr], u[i] <= V[irreversible_reactions_id[i]])
    optimize!(model)

    # Calculating dual variables:

    dualVar = dual.(c1)

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
    return blocked_index, dualVar
end

end
