#-------------------------------------------------------------------------------------------
#=
    Purpose:    Identifying blocked reactions in metabolic networks using linear programming(1LP) and Gaussian elimination
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology
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

The function takes in several arguments, including a ModelObject of type MyModel,
and assigns values to its fields based on the other arguments passed in.

# INPUTS

-'ModelObject':     A newly object of MyModel.
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

The function first exports data from ModelObject and performs some calculations to determine the reversibility of reactions
and the number of blocked reactions. It then creates an optimization model using the GLPK optimizer to identify irreversible blocked reactions,
and uses Gaussian elimination to identify blocked reversible reactions.

# INPUTS

- `ModelObject`:        A newly object of MyModel.

# OPTIONAL INPUTS

- `Tolerance`:          A small number that represents the level of error tolerance.

# OUTPUTS

- `blocked_index`:      IDs of of blocked reactions.
- `dualVar`:            Dual variables of a specific constraint.

# EXAMPLES

- Full input/output example
```julia
julia> blocked_index, dual_var = swiftCC(ModelObject)
```

See also: `MyModel`, myModel_Constructor(), `reversibility()`, `homogenization()`

"""

function swiftCC(ModelObject::MyModel, Tolerance::Float64=1e-6)

    ## Export data from ModelObject

    S = ModelObject.S
    Metabolites = ModelObject.Metabolites
    Reactions = ModelObject.Reactions
    m = ModelObject.m
    n = ModelObject.n
    lb = ModelObject.lb
    ub = ModelObject.ub

    ## Determine the reversibility of a reaction

    irreversible_reactions_id, reversible_reactions_id = reversibility(lb)

    ## Determine the number of irreversible and reversible reactions

    n_irr = length(irreversible_reactions_id)
    n_rev = length(reversible_reactions_id)

    ## Calculate the dimensions of the S matrix

    row_num, col_num = size(S)

    ## Find irreversible blocked reactions in 1 LP problem

    # Set the lower and upper bounds of the "u" variables to 0 and 1, respectively:
    lb_u = zeros(n_irr)
    ub_u = ones(n_irr)

    # Create a new optimization model using the GLPK optimizer:
    model = Model(GLPK.Optimizer)

    # Define variables V and u with their lower and upper bounds:
    @variable(model, lb[i] <= V[i = 1:n] <= ub[i])
    @variable(model, lb_u[i] <= u[i = 1:n_irr] <= ub_u[i])

    # Add a constraint that ensures the mass balance of the system:
    @constraint(model, c1, S * V .== 0)

    # Define the objective function to be the sum of all "u" variables:
    objective_function = ones(n_irr)'*u

    # Set the optimization objective to maximize the objective function:
    @objective(model, Max, objective_function)

    # Add a constraint for each irreversible reaction that sets its "u" variable to be less than or equal to its corresponding V variable:
    @constraint(model, [i in 1:n_irr], u[i] <= V[irreversible_reactions_id[i]])

    # Solve the optimization problem:
    optimize!(model)

    # Get the dual variables for the stoichiometric constraint:
    dualVar = dual.(c1)

    # Initialize an empty list to store the IDs of blocked irreversible reactions:
    irr_blocked_reactions = []

    # Iterate over the irreversible reactions and check if the corresponding u variable is close to 0, If it is, the reaction is considered blocked and its ID is added to the list:
    for i in range(1,n_irr; step=1)
        if isapprox(value(u[i]), 0.0, atol = Tolerance)
            append!(irr_blocked_reactions, irreversible_reactions_id[i])
        end
    end

    # Determine the number of irreversible blocked reactions:
    irr_blocked_num = length(irr_blocked_reactions)

    ## Find reversible blocked reactions

    # Creating S_transpose:
    S_transpose = S'

    # Determining the dimensions of the S_transpose matrix:
    row_num_trans, col_num_trans = size(S_transpose)

    # Removing irreversibly blocked reactions from the Stoichiometric Matrix:
    S_transpose_noIrrBlocked = S_transpose[setdiff(1:end, irr_blocked_reactions), :]

    # Creating the I_reversible Matrix:
    I_reversible = zeros(n, n_rev)

    # Populating the I_reversible Matrix with 1 in the rows corresponding to reversible reactions:
    rev_id = 1
    for col in eachcol(I_reversible)
        col[reversible_reactions_id[rev_id]] = 1.0
        rev_id += 1
    end

    # Removing irreversibly blocked reactions from the I_reversible Matrix:
    I_reversible = I_reversible[setdiff(1:end, irr_blocked_reactions), :]

    # Determining the dimensions of the S_transpose_noIrrBlocked and I_reversible matrices:
    S_trans_row, S_trans_col = size(S_transpose_noIrrBlocked)
    I_row, I_col = size(I_reversible)

    # Solving the system of equations using Gaussian elimination to identify blocked reversible reactions:
    X = S_transpose_noIrrBlocked \ I_reversible
    Sol = (S_transpose_noIrrBlocked*X) - I_reversible

    # Determining the dimensions of the Sol matrix:
    row_sol, col_sol = size(Sol)

    # Finding columns in Sol that correspond to blocked reversible reactions:
    c = 0
    rev_blocked_reactions_col = []
    for col in eachcol(Sol)
        c += 1
        if isapprox(norm(col), 0, atol = Tolerance)
            append!(rev_blocked_reactions_col, c)
        end
    end

    rev_blocked_reactions = []
    for i in rev_blocked_reactions_col
        append!(rev_blocked_reactions, reversible_reactions_id[i])
    end

    ## Union reversbile blocked reactions and irreversible blocked reactions

    blocked_index = []
    blocked_index = union(rev_blocked_reactions, irr_blocked_reactions)

    ## Returning a list consist of the Ids of the blocked reactions

    blocked_index = sort(blocked_index)
    return blocked_index, dualVar
end

end
