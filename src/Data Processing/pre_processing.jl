#-------------------------------------------------------------------------------------------
#=
    Purpose:    A Toolkit of Functions for Preprocessing Metabolic Network Data
    Authors:    Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology
    Date:       April 2022
=#
#-------------------------------------------------------------------------------------------

module pre_processing

export dataOfModel, getM, getTolerance, reversibility, check_duplicate_reaction, homogenization, distributedReversibility_Correction

using GLPK, JuMP, COBREXA, SparseArrays, Distributed, SharedArrays

"""
    dataOfModel(myModel)

The function extracts various data from a given metabolic network model,
represented by a StandardModel object, and returns it as a tuple.

# INPUTS

- `myModel`:        A StandardModel that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

-

# OUTPUTS

- `S`:              A `m` x `n` matrix representing the stoichiometric coefficients of each metabolite in each reaction.
- `Metabolites`:    A list of all metabolites in the metabolic network.
- `Reactions`:      A list of all reactions in the metabolic network.
- `Genes`:          A list of all genes associated with reactions in the metabolic network.
- `m`:              The number of rows in the stoichiometric matrix.
- `n`:              The number of columns in the stoichiometric matrix.
- `lb`:             A `n` x `1` vector representing the lower bounds of each reaction.
- `ub`:             A `n` x `1` vector representing the upper bounds of each reaction.

# EXAMPLES

- Full input/output example
```julia
julia> S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel)
```

"""

function dataOfModel(myModel::StandardModel)

    ## Extracting Data

    S = stoichiometry(myModel) # Stoichiometric matrix
    Metabolites = metabolites(myModel) # Array of metabolite IDs
    Reactions = reactions(myModel) # Array of reaction IDs
    Genes = genes(myModel) # Array of gene IDs
    m = length(metabolites(myModel)) # Number of metabolites
    n = length(reactions(myModel)) # Number of reactions
    lb = lower_bounds(myModel) # Array of lower bounds for each reaction
    ub = upper_bounds(myModel) # Array of upper bounds for each reaction

    ## Sorting Reactions

    p = sortperm(Reactions) # Sort the reaction IDs in alphabetical order
    Reactions = Reactions[p] # Reorder the reaction IDs according to the sorted indices
    lb = lb[p] # Reorder the lower bounds according to the sorted indices
    ub = ub[p] # Reorder the upper bounds according to the sorted indices
    S = S[:,p] # Reorder the columns of the stoichiometric matrix according to the sorted indices

    # Return the extracted data as a tuple:
    return S, Metabolites, Reactions, Genes, m, n, lb, ub
end

#-------------------------------------------------------------------------------------------

"""
    getM()

The function reads the first line of a configuration file named "ConfigFile.txt" located in a subdirectory called "config".
The value on that line is parsed as a floating-point number and returned as the output of the function,
which is typically used to set the upper bound of a variable representing an infinite boundary in a metabolic network model.
If the file cannot be opened or the value on the first line is not a valid floating-point number, the function returns nothing.

# INPUTS

-

# OPTIONAL INPUTS

-

# OUTPUTS

- `M`:             A large number.

# EXAMPLES

- Full input/output example
```julia
julia> M = getM()
```

"""

function getM()

    ## Attempt to open the file "ConfigFile.txt" for reading

    try
        f = open("../config/ConfigFile.txt", "r")

        # Read the first line of the file:
        line = readline(f)

        # Convert the line to a float and close the file:
        M = parse(Float64, line)
        close(f)

        # Return M:
        return M
    catch e
        # If an error occurs, print an error message and return nothing:
        println("Error: Could not open file \"ConfigFile.txt\" for reading.")
        return nothing
    end
end

#-------------------------------------------------------------------------------------------

"""
    getTolerance()

The function reads a value from a configuration file and returns it as a floating-point number to set
the error tolerance level for subsequent computations. If an error occurs while opening or reading the file,
the function prints an error message and returns nothing.

# INPUTS

-

# OPTIONAL INPUTS

-

# OUTPUTS

- `Tolerance`:              A small number that represents the level of error tolerance.

# EXAMPLES

- Full input/output example
```julia
julia> Tolerance = getTolerance()
```

"""

function getTolerance()

    ## Attempt to open the file "ConfigFile.txt" for reading

    try
        f = open("../config/ConfigFile.txt", "r")

        # Read the second line of the file:
        readline(f)
        line = readline(f)

        # Convert the line to a float and close the file:
        Tolerance = parse(Float64, line)
        close(f)

        # Return Tolerance:
        return Tolerance
    catch e
        # If an error occurs, print an error message and return nothing:
        println("Error: Could not open file \"ConfigFile.txt\" for reading.")
        return nothing
    end
end

#-------------------------------------------------------------------------------------------

"""
    reversibility(lb)

The function determines the reversibility of reactions in a metabolic network,
based on their lower bounds. It takes an array of lower bounds as input, and returns two arrays:
one containing the IDs of irreversible reactions, and the other containing the IDs of reversible reactions

# INPUTS

- `lb`:                             Array of lower bounds.

# OPTIONAL INPUTS

-

# OUTPUTS

- `irreversible_reactions_id`:      IDs of irreversible reactions.
- `reversible_reactions_id`:        IDs of reversible reactions.

# EXAMPLES

- Full input/output example
```julia
julia> irreversible_reactions_id, reversible_reactions_id = reversibility(lb)
```

See also: `dataOfModel()`

"""

function reversibility(lb::Array{Float64,1})

    ## Get the length of the "lb" array

    n = length(lb)

    ## Create empty arrays to hold the IDs of irreversible and reversible reactions

    irreversible_reactions_id = Array{Int64}([])
    reversible_reactions_id = Array{Int64}([])

    ## Loop through each reaction in the "lb" array

    for i in 1:n
        # If the lower bound of the reaction is greater than or equal to zero, add the reaction ID to the irreversible reactions array:
        if lb[i] >= 0
            append!(irreversible_reactions_id, i)
        # Otherwise, add the reaction ID to the reversible reactions array:
        else
            append!(reversible_reactions_id, i)
        end
    end
    return irreversible_reactions_id, reversible_reactions_id
end

#-------------------------------------------------------------------------------------------

"""
    check_duplicate_reaction(Reactions)

The function takes an array of strings, representing reaction IDs, and checks if there are any duplicates.

# INPUTS

- `Reactions`:      A list of all reactions in the metabolic network.

# OPTIONAL INPUTS

-

# OUTPUTS

- `true`:           There is a duplicate reactions.
- `false`:          There are no duplicate reactions.

# EXAMPLES

- Full input/output example
```julia
julia> check_duplicate = check_duplicate_reaction(Reactions)
```

See also: `dataOfModel()`

"""

function check_duplicate_reaction(Reactions::Array{String,1})

    # Get the length of the Reactions array:
    n = length(Reactions)
    # Get unique elements of Reactions and modify Reactions to only include unique elements:
    unique_reactions = unique!(Reactions)
    n_unique = length(unique_reactions)

    # Check if the length of Reactions is equal to the length of unique_reactions
    # If they are equal, then there are no duplicate reactions, so return false:
    if n == n_unique
        return false
    else
        return true
    end

end

#-------------------------------------------------------------------------------------------

"""
    homogenization(lb,ub)

The function homogenizes the lower and upper bounds of a set of reactions in a metabolic network model, represented as arrays lb and ub.

# INPUTS

- `lb`:             A `n` x `1` vector representing the lower bounds of each reaction.
- `ub`:             A `n` x `1` vector representing the upper bounds of each reaction.

# OPTIONAL INPUTS

-

# OUTPUTS

- `lb`:             A homogenous `n` x `1` vector representing the lower bounds of each reaction.
- `ub`:             A homogenous `n` x `1` vector representing the upper bounds of each reaction.

# EXAMPLES

- Full input/output example
```julia
julia> lb,ub = homogenization(lb,ub)
```

See also: `dataOfModel()`, 'getM()'

"""

function homogenization(lb::Array{Float64,1}, ub::Array{Float64,1})
    n = length(lb)
    # Set a large number for M:
    M = 1000000.0

    # If the lower bound is greater than zero, set it to zero:
    lb[lb .> 0] .= 0
    # If the upper bound is greater than zero, set it to the constant M:
    ub[ub .> 0] .= M
    # If the lower bound is less than zero, set it to the negative of the constant M:
    lb[lb .< 0] .= -M
    # If the upper bound is less than zero, set it to zero:
    ub[ub .< 0] .= 0

    return lb,ub
end

#-------------------------------------------------------------------------------------------

"""
    reversibility_checking(S, lb, ub, reversible_reactions_id)

A function that detects reversible reactions that are blocked in only one direction.

# INPUTS

- `S`:                           Stoichiometric matrix.
- `lb`:                          LowerBound Of Reactions.
- `ub`:                          UpperBound of Reactions.
- `reversible_reactions_id`:     Reversible reaction IDs.

# OPTIONAL INPUTS

-

# OUTPUTS

- `rev_blocked_fwd`:             Reversible reactions that are blocked in forward direction.
- `rev_blocked_back`:            Reversible reactions that are blocked in backward direction.

# EXAMPLES

- Full input/output example
```julia
julia> rev_blocked_fwd, rev_blocked_back = reversibility_checking(S, lb, ub, reversible_reactions_id)
```

See also: `dataOfModel()`, `reversibility()`, 'getTolerance()'

"""

function reversibility_checking(S::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}, lb::Array{Float64,1}, ub::Array{Float64,1}, reversible_reactions_id::Vector{Int64})

    # Define the number of variables in the model:
    n = length(lb)

    # Set the tolerance value:
    Tolerance = getTolerance()

    # Create a GLPK model object:
    model = Model(GLPK.Optimizer)

    # Define variables V:
    @variable(model, lb[i] <= V[i = 1:n] <= ub[i])

    # Add the stoichiometric constraints to the model:
    @constraint(model, S * V .== 0)

    # Initialize empty arrays to store blocked reactions:
    rev_blocked_fwd = Array{Int64}([])
    rev_blocked_back = Array{Int64}([])

    # Iterate over all reversible reactions in the model

    for j in reversible_reactions_id

        # Set the objective function to maximize the flux through reaction j in the forward direction:
        @objective(model, Max, V[j])

        # Add a constraint that limits the flux through reaction j in the forward direction to be less than or equal to 1:
        @constraint(model, c1, V[j] <= 1)

        # Optimize the model and retrieve the objective value:
        optimize!(model)
        opt_fwd = objective_value(model)

        # If the objective value is approximately 0, the reaction is considered to be blocked in the forward direction:
        if isapprox(opt_fwd, 0, atol=Tolerance)
            append!(rev_blocked_fwd, j)
        end

        # Delete the constraint and unregister it from the model:
        delete(model, c1)
        unregister(model, :c1)

        # Set the objective function to minimize the flux through reaction j in the backward direction:
        @objective(model, Min, V[j])

        # Add a constraint that limits the flux through reaction j in the backward direction to be greater than or equal to -1:
        @constraint(model, c2, V[j] >= -1)

        # Optimize the model and retrieve the objective value:
        optimize!(model)
        opt_back = objective_value(model)

        # If the objective value is approximately 0, the reaction is considered to be blocked in the backward direction:
        if isapprox(opt_back, 0, atol=Tolerance)
            append!(rev_blocked_back, j)
        end

        # Delete the constraint and unregister it from the model:
        delete(model, c2)
        unregister(model, :c2)
    end

    return rev_blocked_fwd, rev_blocked_back
end

#-------------------------------------------------------------------------------------------

"""
    reversibility_correction(S, lb, ub, irreversible_reactions_id, reversible_reactions_id, rev_blocked_fwd, rev_blocked_back)

A function that modifies 3 sets:

    1) Remove rev_blocked_fwd and rev_blocked_back from reversible reactions list.
    2) Add rev_blocked_fwd and rev_blocked_back to irreversible reactions list.
    3) modify S, lb and ub.

# INPUTS

- `S`:                                    Stoichiometric matrix.
- `lb`:                                   LowerBound Of Reactions.
- `ub`:                                   UpperBound of Reactions.
- `irreversible_reactions_id`:            Irreversible reaction IDs.
- `reversible_reactions_id`:              Reversible reaction IDs.
- `rev_blocked_fwd`:                      Reversible reactions that are blocked in forward direction.
- `rev_blocked_back`:                     Reversible reactions that are blocked in backward direction.

# OPTIONAL INPUTS

-

# OUTPUTS

- `S`:                                    Stoichiometric matrix.
- `lb`:                                   LowerBound Of Reactions.
- `ub`:                                   UpperBound of Reactions.
- `irreversible_reactions_id`:            Irreversible reaction IDs.
- `reversible_reactions_id`:              Corrected reversible reaction IDs.

# EXAMPLES

- Full input/output example
```julia
julia> S, lb, ub, irreversible_reactions_id, reversible_reactions_id = reversibility_correction(S, lb, ub, irreversible_reactions_id, reversible_reactions_id, rev_blocked_fwd, rev_blocked_back)
```

See also: `dataOfModel()`, `homogenization()`, `reversibility()`, reversibility_checking()

"""

function reversibility_correction(S::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}, lb::Array{Float64,1}, ub::Array{Float64,1}, irreversible_reactions_id::Vector{Int64}, reversible_reactions_id::Vector{Int64}, rev_blocked_fwd::Vector{Int64}, rev_blocked_back::Vector{Int64})

    corrected_reversible_reactions_id = Array{Int64}([])

    ## Forward

    for i in rev_blocked_fwd
        # Modify lower and upper bounds:
        ub[i] = lb[i] * -1
        lb[i] = 0.0
        # Add to irreversible reactions list:
        append!(irreversible_reactions_id, i)
        # Modify Stoichiometric Matrix:
        S[:, i] .= S[:, i] * -1
    end

    ## Backward

    for i in rev_blocked_back
        # Modify lower bounds:
        lb[i] = 0.0
        # Add to irreversible reactions list:
        append!(irreversible_reactions_id, i)
    end

    # Remove blocked reversible reactions from the reversible reactions list

    set_reversible_reactions_id = Set(reversible_reactions_id)
    set_rev_blocked_fwd = Set(rev_blocked_fwd)
    set_rev_blocked_back = Set(rev_blocked_back)
    set_rev_blocked_onedirection = union(set_rev_blocked_fwd, set_rev_blocked_back)
    set_reversible_reactions_id = setdiff(set_reversible_reactions_id, set_rev_blocked_onedirection)

    # Add remaining reversible reactions to the corrected reversible reactions list

    for i in set_reversible_reactions_id
        append!(corrected_reversible_reactions_id, i)
    end

    return S, lb, ub, irreversible_reactions_id, corrected_reversible_reactions_id
end

end
