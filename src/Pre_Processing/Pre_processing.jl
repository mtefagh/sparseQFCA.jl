#-----------------------------------------------------------------------------------------------------------------------------------------------------
#=
    Purpose:    A Toolkit of Functions for Preprocessing Metabolic Network Data
    Authors:    Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology
    Date:       April 2022
=#
#-----------------------------------------------------------------------------------------------------------------------------------------------------

module Pre_processing

export dataOfModel, getM, getTolerance, reversibility, check_duplicate_reactions, homogenization, remove_zeroRows, Model_Correction, model_Correction_Constructor, distributedReversibility_Correction

using JuMP, COBREXA, SparseArrays, Distributed, SharedArrays, Distributed, Clarabel

include("Solve.jl")

using .Solve

import CDDLib

import AbstractFBCModels as A

"""
    dataOfModel(model)

The function extracts various data from a given metabolic network model, represented by a CanonicalModel object, and returns it as a tuple.

# INPUTS

- `model`:          A CanonicalModel that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

- `printLevel`:     Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `S`:              A `m` x `n` matrix representing the stoichiometric coefficients of each metabolite in each reaction.
- `Metabolites`:    A list of all metabolites in the metabolic network.
- `Reactions`:      A list of all reactions in the metabolic network.
- `Genes`:          A list of all genes associated with reactions in the metabolic network.
- `m`:              The number of rows in the stoichiometric matrix.
- `n`:              The number of columns in the stoichiometric matrix.
- `n_genes`:        The number of columns in the stoichiometric matrix.
- `lb`:             A `n` x `1` vector representing the lower bounds of each reaction.
- `ub`:             A `n` x `1` vector representing the upper bounds of each reaction.
- `c_vector`:       A `n` x `1` vector representing the Objective coefficients vector.

# EXAMPLES

- Full input/output example
```julia
julia> S, Metabolites, Reactions, Genes, m, n, n_genes, lb, ub, c_vector = dataOfModel(model)
```

See also: `COBREXA.load_model(JSONFBCModel, "", A.CanonicalModel.Model)`

"""

function dataOfModel(model, printLevel::Int=1)

    ## Extracting Data

    # Stoichiometric matrix:
    S = A.stoichiometry(model)

    # Array of metabolite IDs:
    Metabolites = A.metabolites(model)

    # Array of reaction IDs:
    Reactions = A.reactions(model)

    # Array of gene IDs:
    Genes = A.genes(model)

    # Number of metabolites in the model:
    m = A.n_metabolites(model)

    # Number of reactions in the model:
    n = A.n_reactions(model)

    # Number of genes in the model:
    n_genes = A.n_genes(model)

    # Array of lower and upper bounds for each reaction (constraints on flux through each reaction):
    lb, ub = A.bounds(model)

    # Objective coefficients vector:
    c_vector = A.objective(model)

    ## Sorting Reactions (Commented out)

    #=
    # Sort reaction IDs in alphabetical order and reorder the reactions and associated data accordingly
    p = sortperm(Reactions)  # Get the sorted indices of the reactions
    Reactions = Reactions[p]  # Reorder reactions by sorted indices
    lb = lb[p]  # Reorder lower bounds accordingly
    ub = ub[p]  # Reorder upper bounds accordingly
    S = S[:,p]  # Reorder columns of the stoichiometric matrix accordingly
    =#

    ## Print out results if requested
    if printLevel > 0
        # Print a header for the output in cyan color
        printstyled("Metabolic Network:\n"; color=:cyan)

        # Print the number of metabolites, reactions, and genes in the model
        println("Number of Metabolites      : $m")
        println("Number of Reactions        : $n")
        println("Number of Genes            : $n_genes")

        # Print the dimensions of the stoichiometric matrix (rows = metabolites, columns = reactions)
        println("Stoichiometric matrix      : $m x $n")
    end

    # Return the extracted data as a tuple
    return S, Metabolites, Reactions, Genes, m, n, n_genes, lb, ub, c_vector
end

#-------------------------------------------------------------------------------------------

"""
    getM()

The function reads the first line of a configuration file named "ConfigFile.txt" located in a subdirectory called "config".
The value on that line is parsed as a floating-point number and returned as the output of the function, which is typically
used to set the upper bound of a variable representing an infinite boundary in a metabolic network model. If the file cannot
be opened or the value on the first line is not a valid floating-point number, the function returns nothing.

# INPUTS

-

# OPTIONAL INPUTS

- `printLevel`:    Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `M`:             A large number.

# EXAMPLES

- Full input/output example
```julia
julia> M = getM()
```

"""

function getM(printLevel::Int=1)

    ## Attempt to open the file "ConfigFile.txt" for reading

    try
        f = open("../config/ConfigFile.txt", "r")

        # Read the first line of the file:
        line = readline(f)

        # Convert the line to a float and close the file:
        M = parse(Float64, line)
        close(f)

        ## Print out results if requested

        if printLevel > 0
            printstyled("M = $M\n"; color=:magenta)
        end

        # Return M:
        return M
    catch e
        # If an error occurs, print an error message and return nothing:
        printstyled("Error: Could not open file \"ConfigFile.txt\" for reading.\n"; color=:red)
        return nothing
    end
end

#-------------------------------------------------------------------------------------------

"""
    getTolerance()

The function reads a value from a configuration file and returns it as a floating-point number to set the
error tolerance level for subsequent computations. If an error occurs while opening or reading the file,the
function prints an error message and returns nothing.

# INPUTS

-

# OPTIONAL INPUTS

- `printLevel`:             Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `Tolerance`:              A small number that represents the level of error tolerance.

# EXAMPLES

- Full input/output example
```julia
julia> Tolerance = getTolerance()
```

"""

function getTolerance(printLevel::Int=0)

    ## Attempt to open the file "ConfigFile.txt" for reading

    try
        f = open("../config/ConfigFile.txt", "r")

        # Read the second line of the file:
        readline(f)
        line = readline(f)

        # Convert the line to a float and close the file:
        Tolerance = parse(Float64, line)
        close(f)

        ## Print out results if requested

        if printLevel > 0
            printstyled("Tolerance = $Tolerance\n"; color=:magenta)
        end

        # Return Tolerance:
        return Tolerance
    catch e
        # If an error occurs, print an error message and return nothing:
        printstyled("Error: Could not open file \"ConfigFile.txt\" for reading.\n"; color=:red)
        return nothing
    end
end

#-------------------------------------------------------------------------------------------

"""
    reversibility(lb, Reaction_Ids)

The function determines the reversibility of reactions in a metabolic network,
based on their lower bounds. It takes an array of lower bounds as input, and returns two arrays:
one containing the IDs of irreversible reactions, and the other containing the IDs of reversible reactions

# INPUTS

- `lb`:                             Array of lower bounds.
- `Reaction_Ids`:                   Array of reaction IDs.

# OPTIONAL INPUTS

- `printLevel`  :                   Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `irreversible_reactions_id`:      IDs of irreversible reactions.
- `reversible_reactions_id`  :      IDs of reversible reactions.

# EXAMPLES

- Full input/output example
```julia
julia> irreversible_reactions_id, reversible_reactions_id = reversibility(lb, Reaction_Ids)
```

See also: `dataOfModel()`

"""

function reversibility(lb::Array{Float64,1}, Reaction_Ids::Vector{Int64}, printLevel::Int=1)

    # Get the length of the "lb" array:
    n = length(lb)

    # Create empty arrays to hold the IDs of irreversible and reversible reactions:
    irreversible_reactions_id = Array{Int64}([])
    reversible_reactions_id = Array{Int64}([])

    ## Loop through each reaction in the "lb" array

    for i = 1:n
        # If the lower bound of the reaction is greater than or equal to zero, add the reaction ID to the irreversible reactions array:
        if lb[i] >= 0
            append!(irreversible_reactions_id, Reaction_Ids[i])
        # Otherwise, add the reaction ID to the reversible reactions array:
        else
            append!(reversible_reactions_id, Reaction_Ids[i])
        end
    end

    # Get the number of irreversible and reversible reactions:
    n_irr = length(irreversible_reactions_id)
    n_rev = length(reversible_reactions_id)

    ## Print out results if requested

    if printLevel > 0
        printstyled("Reversibility Checking:\n"; color=:cyan)
        println("Number of irreversible reactions : $n_irr ")
        println("Number of reversible   reactions : $n_rev ")
    end

    # Return the IDs of the irreversible and reversible reactions:
    return irreversible_reactions_id, reversible_reactions_id
end

#-------------------------------------------------------------------------------------------

"""
    check_duplicate_reaction(Reactions)

The function takes an array of strings, representing reaction IDs, and checks if there are any duplicates.

# INPUTS

- `Reactions`:      A list of all reactions in the metabolic network.

# OPTIONAL INPUTS

- `printLevel`:     Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `true`:           There is a duplicate reaction.
- `false`:          There is no duplicate reaction.

# EXAMPLES

- Full input/output example
```julia
julia> check_duplicate_reactions(Reactions)
```

See also: `dataOfModel()`

"""

function check_duplicate_reactions(Reactions::Array{String,1}, printLevel::Int=1)

    # Get the length of the Reactions array:
    n = length(Reactions)
    # Get unique elements of Reactions and modify Reactions to only include unique elements:
    unique_reactions = unique!(Reactions)
    n_unique = length(unique_reactions)

    # Check if the length of Reactions is equal to the length of unique_reactions
    # If they are equal, then there are no duplicate reactions, so return false:
    if n == n_unique
        printstyled("There is no duplicate reaction.\n"; color=:green)
        return false
    else
        printstyled("There is a duplicate reaction.\n"; color=:red)
        return true
    end
end

#-------------------------------------------------------------------------------------------

"""
    homogenization(lb, ub)

The function homogenizes the lower and upper bounds of a set of reactions in a metabolic network model, represented as arrays lb and ub.

# INPUTS

- `lb`:             A `n` x `1` vector representing the lower bounds of each reaction.
- `ub`:             A `n` x `1` vector representing the upper bounds of each reaction.

# OPTIONAL INPUTS

- `printLevel`:     Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `lb`:             A homogenous `n` x `1` vector representing the lower bounds of each reaction.
- `ub`:             A homogenous `n` x `1` vector representing the upper bounds of each reaction.

# EXAMPLES

- Full input/output example
```julia
julia> lb, ub = homogenization(lb, ub)
```

See also: `dataOfModel()`, 'getM()'

"""

function homogenization(lb::Array{Float64,1}, ub::Array{Float64,1}, printLevel::Int=1)

    # Get the length of the lb array:
    n = length(lb)

    if printLevel > 0
        printstyled("Homogenization:\n"; color=:cyan)
    end

    # Set a large number for M:
    M = getM(printLevel)

    # If the lower bound is greater than zero, set it to zero:
    lb[lb .>= 0] .= 0
    # If the upper bound is greater than zero, set it to the constant M:
    ub[ub .> 0] .= M
    # If the lower bound is less than zero, set it to the negative of the constant M:
    lb[lb .< 0] .= -M
    # If the upper bound is less than zero, set it to zero:
    ub[ub .<= 0] .= 0

    # Return the homogenized lower and upper bounds as a tuple:
     return lb, ub
end

#-------------------------------------------------------------------------------------------

"""
    remove_zeroRows(S, Metabolites)

The function returns a modified stoichiometric matrix S and corresponding array of metabolites Metabolites,
with all rows that contain only zeros removed. This results in a more compact representation of the metabolic
network, with only the relevant reactions and metabolites included.

# INPUTS

- `S`:                      A `m` x `n` matrix representing the stoichiometric coefficients of each metabolite in each reaction.
- `Metabolites`:            A list of all metabolites in the metabolic network.

# OPTIONAL INPUTS

-

# OUTPUTS

- `S`:                      The stoichiometry matrx of non-zero interactions within the metabolic network.
- `Metabolites`:            The updated list of all metabolites in the metabolic network.

# EXAMPLES

- Full input/output example
```julia
julia> S, Metabolites, Metabolites_elimination = remove_zeroRows(S, Metabolites)
```

See also: `dataOfModel()`

"""

function remove_zeroRows(S::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}, Metabolites::Array{String,1})

    zero_row = [] # create an empty array to store indices of zero rows
    c = 1 # initialize a counter variable to keep track of row indices
    for row ∈ eachrow(S) # iterate over each row of the matrix S
        if row == zeros(length(row)) # check if the current row is all zeros
            append!(zero_row, c) # if so, append the current row index to zero_row array
        end
        c += 1 # increment the row index counter
    end

    # Remove the zero rows from the matrix S:
    S = S[setdiff(1:end, zero_row), :]

    # get the length of the Metabolites array:
    m = length(Metabolites)

    # Remove the corresponding metabolites from the Metabolites array:
    Metabolites_elimination = Metabolites[zero_row]
    Metabolites = Metabolites[setdiff(range(1, m), zero_row)]

    # Return the updated S matrix, Metabolites array and Metabolites_elimination array:
    return S, Metabolites, Metabolites_elimination
end

#-------------------------------------------------------------------------------------------

"""
    Model_Correction(S, Metabolites, Reactions, Genes, m, n, lb, ub, irreversible_reactions_id, reversible_reactions_id)

A general type for storing a CoreModel which contains the following fields to run distributedReversibility_Correction():

- `S`:                              LHS matrix (m x n)
- `Metabolites`:                    List of metabolic network metabolites.
- `Reactions`:                      List of metabolic network reactions.
- `Genes`:                          List of metabolic network reactions.
- `m`:                              Number of rows of stoichiometric matrix.
- `n`:                              Number of columns of stoichiometric matrix.
- `lb`:                             Lower bound vector (n x 1)
- `ub`:                             Upper bound vector (n x 1)
- `irreversible_reactions_id`:      IDs of irreversible reactions.
- `reversible_reactions_id`  :      IDs of reversible reactions.

"""

mutable struct Model_Correction
    S                               ::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}
    Metabolites                     ::Array{String,1}
    Reactions                       ::Array{String,1}
    Genes                           ::Array{String,1}
    m                               ::Int
    n                               ::Int
    lb                              ::Array{Float64,1}
    ub                              ::Array{Float64,1}
    irreversible_reactions_id       ::Array{Int64}
    reversible_reactions_id         ::Array{Int64}
end

#-------------------------------------------------------------------------------------------

"""
    model_Correction_Constructor(ModelObject_Correction, S, Metabolites, Reactions, Genes, m, n, lb, ub, irreversible_reactions_id, reversible_reactions_id)

The function takes in several arguments, including a ModelObject of type Model_Correction,
and assigns values to its fields based on the other arguments passed in.

# INPUTS

-`ModelObject_Correction`:          A newly object of Model_Correction.
- `S`:                              LHS matrix (m x n)
- `Metabolites`:                    List of metabolic network metabolites.
- `Reactions`:                      List of metabolic network reactions.
- `Genes`:                          List of metabolic network reactions.
- `m`:                              Number of rows of stoichiometric matrix.
- `n`:                              Number of columns of stoichiometric matrix.
- `lb`:                             Lower bound vector (n x 1)
- `ub`:                             Upper bound vector (n x 1)
- `irreversible_reactions_id`:      IDs of irreversible reactions.
- `reversible_reactions_id`  :      IDs of reversible reactions.

"""

function model_Correction_Constructor(ModelObject_Correction::Model_Correction, S::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}, Metabolites::Array{String,1}, Reactions::Array{String,1},
                                      Genes::Array{String,1},m::Int, n::Int, lb::Array{Float64,1}, ub::Array{Float64,1}, irreversible_reactions_id::Array{Int64}, reversible_reactions_id::Array{Int64})
     ModelObject_Correction.S = S
     ModelObject_Correction.Metabolites = Metabolites
     ModelObject_Correction.Reactions = Reactions
     ModelObject_Correction.Genes = Genes
     ModelObject_Correction.m = m
     ModelObject_Correction.n = n
     ModelObject_Correction.lb = lb
     ModelObject_Correction.ub = ub
     ModelObject_Correction.irreversible_reactions_id = irreversible_reactions_id
     ModelObject_Correction.reversible_reactions_id = reversible_reactions_id
end

#-------------------------------------------------------------------------------------------

"""
    distributedReversibility_Correction(ModelObject_Correction, blocked_index_rev)

The function first distributes the list of reversible reactions among worker processes, and for each reversible reaction,
solves two linear programming problems to maximize the flux through the reaction in the forward direction and minimize the flux
through the reaction in the backward direction. The function then identifies which reactions are blocked in each direction,
modifies the lower and upper bounds accordingly, and appends the blocked reactions to the list of irreversible reactions.
Finally, the function removes the blocked reactions from the list of reversible reactions.

# INPUTS

-`ModelObject_Correction`:                  A newly object of Model_Correction.
- `blocked_index`:                          IDs of reversible blocked reactions.

# OPTIONAL INPUTS

- `SolverName`:                             Name of the solver(default: HiGHS).
- `OctuplePrecision`:                       A flag(default: false) indicating whether octuple precision should be used when solving linear programs.
- `printLevel`:                             Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `S`:                                      Stoichiometric matrix with corrected reversible reactions.
- `lb`:                                     A corrected `n` x `1` vector representing the lower bounds of each reaction.
- `ub`:                                     A corrected `n` x `1` vector representing the upper bounds of each reaction.
- `irreversible_reactions_id`:              IDs of corrected irreversible reactions.
- `corrected_reversible_reactions_id`:      IDs of corrected reversible reactions.

# EXAMPLES

- Full input/output example
```julia
julia> S, lb, ub, irreversible_reactions_id, reversible_reactions_id = distributedReversibility_Correction(ModelObject_Correction, blocked_index_rev)
```

See also: `dataOfModel()`, `reversibility()`, `homogenization()`, 'getTolerance()'

"""

function distributedReversibility_Correction(ModelObject_Correction::Model_Correction, blocked_index_rev::Vector{Int64}, SolverName::String="HiGHS", OctuplePrecision::Bool=false, printLevel::Int=1)

    ## Extract relevant information from the ModelObject_Correction

    # Extracting the stoichiometric matrix (S) from the ModelObject_Correction:
    S = ModelObject_Correction.S

    # Extracting the array of metabolite IDs from the ModelObject_Correction:
    Metabolites = ModelObject_Correction.Metabolites

    # Extracting the array of reaction IDs from the ModelObject_Correction:
    Reactions = ModelObject_Correction.Reactions

    # Extracting the array of gene IDs from the ModelObject_Correction:
    Genes = ModelObject_Correction.Genes

    # Extracting the number of metabolites (m) from the ModelObject_Correction:
    m = ModelObject_Correction.m

    # Extracting the number of reactions (n) from the ModelObject_Correction:
    n = ModelObject_Correction.n

    # Extracting the lower bounds for reactions from the ModelObject_Correction:
    lb = ModelObject_Correction.lb

    # Extracting the upper bounds for reactions from the ModelObject_Correction:
    ub = ModelObject_Correction.ub

    # Extracting the IDs of irreversible reactions from the ModelObject_QFCA:
    irreversible_reactions_id = ModelObject_Correction.irreversible_reactions_id

    # Extracting the IDs of reversible reactions from the ModelObject_QFCA:
    reversible_reactions_id = ModelObject_Correction.reversible_reactions_id

    lb_homo = copy(lb)
    ub_homo = copy(ub)

    lb_homo, ub_homo = homogenization(lb_homo, ub_homo, 0)

    # Define the number of variables in the model:
    n = length(lb)

    # Set the tolerance value:
    Tolerance = getTolerance()

    # Initialize empty arrays to store the IDs of blocked reversible reactions in the forward and backward directions:
    rev_blocked_fwd = Array{Int64}([])
    rev_blocked_back = Array{Int64}([])

    # Initialize a shared array to represent whether each reversible reaction is blocked in the forward or backward direction
    # The initial value of false indicates that no reactions are blocked at the beginning:
    Correction = SharedArray{Int,2}((n, n), init = false)

    # Create a local model object for each worker process:
    if OctuplePrecision
        # Define a local model using GenericModel from Clarabel.jl:
        local_model = GenericModel{BigFloat}(Clarabel.Optimizer{BigFloat})
        # Set verbose attribute to false (disable verbose output):
        set_attribute(local_model, "verbose", false)
        # Set absolute tolerance for gap convergence to 1e-32:
        set_attribute(local_model, "tol_gap_abs", 1e-32)
        # Set relative tolerance for gap convergence to 1e-32:
        set_attribute(local_model, "tol_gap_rel", 1e-32)
    else
        local_model, solver = changeSparseQFCASolver(SolverName)
    end

    # Define variables V:
    @variable(local_model, lb_homo[i] <= V[i = 1:n] <= ub_homo[i])

    # Add the stoichiometric constraints to the model:
    @constraint(local_model, S * V .== 0)

    # Iterate over all reversible reactions in the model:
    @sync @distributed for i ∈ reversible_reactions_id

        if i ∈ blocked_index_rev
            # The reaction is considered to be blocked:
            Correction[i,i] = +2.0
        else
            # Set the objective function to maximize the flux through reaction i in the forward direction:
            @objective(local_model, Max, V[i])

            # Add a constraint that limits the flux through reaction i in the forward direction to be less than or equal to 1:
            @constraint(local_model, c_irr, V[i] <= 1)

            # Optimize the model and retrieve the objective value:
            optimize!(local_model)
            opt_fwd = objective_value(local_model)

            # Remove the current constraint from the model:
            delete(local_model, c_irr)

            # Unregister the current constraint from the model:
            unregister(local_model, :c_irr)

            # Set the objective function to minimize the flux through reaction i in the backward direction:
            @objective(local_model, Min, V[i])

            # Add a constraint that limits the flux through reaction i in the backward direction to be greater than or equal to -1:
            @constraint(local_model, c_rev, V[i] >= -1)

            # Optimize the model and retrieve the objective value:
            optimize!(local_model)
            opt_back = objective_value(local_model)

            # Remove the current constraint from the model:
            delete(local_model, c_rev)

            # Unregister the current constraint from the model:
            unregister(local_model, :c_rev)

            # If the objective value is approximately 0, the reaction is considered to be blocked in the forward direction:
            if isapprox(opt_fwd, 0, atol=Tolerance)
                Correction[i,i] = +1.0
            # If the objective value is approximately 0, the reaction is considered to be blocked in the backward direction:
            elseif isapprox(opt_back, 0, atol=Tolerance)
                Correction[i,i] = -1.0
            else
                continue
            end
        end
    end

    for i = 1 : n
        for j = 1 : n
            # If the (i,j)-th element of Correction is +1.0 and i == j,
            # then the i-th reversible reaction is blocked in the forward direction:
            if Correction[i,j] == +1.0 && i == j
                append!(rev_blocked_fwd, i)
            # If the (i,j)-th element of Correction is -1.0 and i == j,
            # then the i-th reversible reaction is blocked in the backward direction:
            elseif Correction[i,j] == -1.0 && i == j
                append!(rev_blocked_back, i)
            # If neither of the above conditions are met, continue to the next iteration of the inner loop:
            else
                continue
            end
        end
    end

    # Initialize an empty array to store the IDs of corrected reversible reactions:
    corrected_reversible_reactions_id = Array{Int64}([])

    ## Forward

    for i ∈ rev_blocked_fwd
        # Modify lower and upper bounds:
        ub[i] = lb[i] * -1
        lb[i] = 0.0
        # Add to irreversible reactions list:
        append!(irreversible_reactions_id, i)
        # Modify Stoichiometric Matrix:
        S[:, i] .= S[:, i] * -1
    end

    ## Backward

    for i ∈ rev_blocked_back
        # Modify lower bounds:
        lb[i] = 0.0
        # Add to irreversible reactions list:
        append!(irreversible_reactions_id, i)
    end

    ## Remove blocked reversible reactions from the reversible reactions list

    # Convert the lists to sets to perform set operations:
    set_reversible_reactions_id = Set(reversible_reactions_id)
    set_rev_blocked_fwd = Set(rev_blocked_fwd)
    set_rev_blocked_back = Set(rev_blocked_back)
    # Find the union of the sets for the reactions blocked in either direction:
    set_rev_blocked_onedirection = union(set_rev_blocked_fwd, set_rev_blocked_back)
    # Remove the blocked reactions from the set of reversible reactions:
    set_reversible_reactions_id = setdiff(set_reversible_reactions_id, set_rev_blocked_onedirection)

    ## Add remaining reversible reactions to the corrected reversible reactions list

    # Convert the remaining reversible reactions set back to a list:
    for i ∈ set_reversible_reactions_id
        append!(corrected_reversible_reactions_id, i)
    end

    ## Print out results if requested

    if printLevel > 0
        printstyled("Distributed Reversibility Correction:\n"; color=:cyan)
        println("Number of Proccess : $(nprocs())")
        println("Number of Workers  : $(nworkers())")
        if OctuplePrecision
            printstyled("Solver = Clarabel \n"; color=:green)
        else
            printstyled("Solver = $SolverName\n"; color=:green)
        end
        printstyled("Tolerance = $Tolerance\n"; color=:magenta)
        println("Number of reversibe blocked in forward  direction : $(length(rev_blocked_fwd))")
        println("Number of reversibe blocked in backward direction : $(length(rev_blocked_back))")
        println("Number of irreversible reactions after Correction : $(length(irreversible_reactions_id))")
        println("Number of reversible   reactions after Correction : $(length(corrected_reversible_reactions_id))")
    end

    # Return the modified model and reaction information as a tuple:
    return S, lb, ub, irreversible_reactions_id, corrected_reversible_reactions_id
end

end

#-----------------------------------------------------------------------------------------------------------------------------------------------------
