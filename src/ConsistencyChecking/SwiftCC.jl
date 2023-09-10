#-----------------------------------------------------------------------------------------------------------------------------------------------------
#=
    Purpose:    Identifying blocked reactions in metabolic networks using Convex Optimization(1LP) and Gaussian Elimination
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology
    Date:       July 2022
=#
#-----------------------------------------------------------------------------------------------------------------------------------------------------

module SwiftCC

export Model_CC, model_CC_Constructor, swiftCC

using GLPK, JuMP, COBREXA, LinearAlgebra, SparseArrays, Distributed

include("../Pre_Processing/Pre_processing.jl")

using .Pre_processing

"""
    Model_CC(S, Metabolites, Reactions, Genes, m, n, lb, ub)

A general type for storing a StandardModel which contains the following fields to run swiftCC():

- `S`:              LHS matrix (m x n)
- `Metabolites`:    List of metabolic network metabolites.
- `Reactions`:      List of metabolic network reactions.
- `Genes`:          List of metabolic network reactions.
- `m`:              Number of rows of stoichiometric matrix.
- `n`:              Number of columns of stoichiometric matrix.
- `lb`:             Lower bound vector (n x 1)
- `ub`:             Upper bound vector (n x 1)

"""

mutable struct Model_CC
    S              ::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}
    Metabolites    ::Array{String,1}
    Reactions      ::Array{String,1}
    Genes          ::Array{String,1}
    m              ::Int
    n              ::Int
    lb             ::Array{Float64,1}
    ub             ::Array{Float64,1}
end

#-------------------------------------------------------------------------------------------

"""
    model_CC_Constructor(ModelObject_CC, S, Metabolites, Reactions, Genes, m, n, lb, ub)

The function takes in several arguments, including a ModelObject of type Model_CC,
and assigns values to its fields based on the other arguments passed in.

# INPUTS

-'ModelObject_CC':     A newly object of Model_CC.
- `S`:                 LHS matrix (m x n)
- `Metabolites`:       List of metabolic network metabolites.
- `Reactions`:         List of metabolic network reactions.
- `Genes`:             List of metabolic network reactions.
- `m`:                 Number of rows of stoichiometric matrix.
- `n`:                 Number of columns of stoichiometric matrix.
- `lb`:                Lower bound vector (n x 1)
- `ub`:                Upper bound vector (n x 1)

"""

function model_CC_Constructor(ModelObject_CC::Model_CC, S::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}, Metabolites::Array{String,1},
                              Reactions::Array{String,1}, Genes::Array{String,1},m::Int, n::Int, lb::Array{Float64,1}, ub::Array{Float64,1})
     ModelObject_CC.S = S
     ModelObject_CC.Metabolites = Metabolites
     ModelObject_CC.Reactions = Reactions
     ModelObject_CC.Genes = Genes
     ModelObject_CC.m = m
     ModelObject_CC.n = n
     ModelObject_CC.lb = lb
     ModelObject_CC.ub = ub
end

#-------------------------------------------------------------------------------------------

"""
    swiftCC(ModelObject_CC)

The function first exports data from ModelObject_CC and performs some calculations to determine the reversibility of reactions
and the number of blocked reactions. It then creates an optimization model using the GLPK optimizer to identify irreversible blocked reactions,
and uses Gaussian elimination to identify blocked reversible reactions.

# INPUTS

- `ModelObject_CC`:        A newly object of Model_CC.

# OPTIONAL INPUTS

- `Tolerance`:             A small number that represents the level of error tolerance.
- `printLevel`:            Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `blocked_index`:         IDs of of blocked reactions.
- `ν`:                     Dual variables of a specific constraint.

# EXAMPLES

- Full input/output example
```julia
julia> blocked_index, ν = swiftCC(ModelObject_CC)
```

See also: `Model_CC`, `model_CC_Constructor()`, `reversibility()`

"""

function swiftCC(ModelObject_CC::Model_CC, Tolerance::Float64=1e-6, printLevel::Int=1)

    ## Extract relevant information from the input model object

    S = ModelObject_CC.S
    Metabolites = ModelObject_CC.Metabolites
    Reactions = ModelObject_CC.Reactions
    m = ModelObject_CC.m
    n = ModelObject_CC.n
    lb = ModelObject_CC.lb
    ub = ModelObject_CC.ub

    ## Identify which reactions are irreversible and which are reversible

    # Create an array of reaction IDs:
    Reaction_Ids = collect(1:n)
    irreversible_reactions_id, reversible_reactions_id = reversibility(lb, Reaction_Ids, 0)
    n_irr = length(irreversible_reactions_id)
    n_rev = length(reversible_reactions_id)

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
    ν = dual.(c1)

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

    # Create S_transpose:
    S_transpose = S'

    # Remove irreversibly blocked reactions from the Stoichiometric Matrix:
    S_transpose_noIrrBlocked = S_transpose[setdiff(1:end, irr_blocked_reactions), :]

    # Create the I_reversible Matrix:
    I_reversible = zeros(n, n_rev)

    # Populate the I_reversible Matrix with 1 in the rows corresponding to reversible reactions:
    rev_id = 1
    for col in eachcol(I_reversible)
        col[reversible_reactions_id[rev_id]] = 1.0
        rev_id += 1
    end

    # Remove irreversibe blocked reactions from the I_reversible Matrix:
    I_reversible = I_reversible[setdiff(1:end, irr_blocked_reactions), :]

    # Determine the dimensions of the S_transpose_noIrrBlocked and I_reversible matrices:
    S_trans_row, S_trans_col = size(S_transpose_noIrrBlocked)
    I_row, I_col = size(I_reversible)

    # Solve the system of equations using Gaussian elimination to identify blocked reversible reactions:
    X = S_transpose_noIrrBlocked \ I_reversible
    Sol = (S_transpose_noIrrBlocked*X) - I_reversible

    # Determine the dimensions of the Sol matrix:
    row_sol, col_sol = size(Sol)

    # Find columns in Sol that correspond to blocked reversible reactions:
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

    ## Print out results if requested

    if printLevel > 0
        printstyled("Consistency_Checking(SwiftCC):\n"; color=:cyan)
        println("Number of Proccess : $(nprocs())")
        println("Number of Workers  : $(nworkers())")
        printstyled("Tolerance = $Tolerance\n"; color=:magenta)
        println("Number of irreversible blocked reactions : $(length(irr_blocked_reactions))")
        println("Number of reversible   blocked reactions : $(length(rev_blocked_reactions))")
        println("Number of blocked reactions              : $(length(blocked_index))")
    end

    return blocked_index, ν
end

end

#-----------------------------------------------------------------------------------------------------------------------------------------------------
