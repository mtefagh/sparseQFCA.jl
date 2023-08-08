#-----------------------------------------------------------------------------------------------------------------------------------------------------
#=
    Purpose:    Identifying Blocked Reactions in a Metabolic Model using Convex Optimization(n_i + 2n_r LP)
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology
    Date:       April 2022
=#
#-----------------------------------------------------------------------------------------------------------------------------------------------------

module TheNaiveApproach

export find_blocked_reactions

using GLPK, JuMP, COBREXA, Distributed

include("../Pre_Processing/Pre_processing.jl")

using .Pre_processing

"""
    find_blocked_reactions(myModel)

The function identifies blocked reactions in a metabolic model. Blocked reactions are reactions that cannot carry any flux,
i.e., the flux through the reaction is zero, due to the stoichiometry of the model and the constraints on the reaction rates.
The function first homogenizes the bounds on the reaction rates, and then identifies blocked irreversible reactions by optimizing
the flux through each reaction while constraining the flux to be less than or equal to 1. If the optimal flux is approximately zero,
the reaction is considered blocked. Then, the function identifies blocked reversible reactions by optimizing the flux through each
reaction in both directions, forward and backward, while constraining the flux to be less than or equal to 1 in the forward direction
and greater than or equal to -1 in the backward direction. If the optimal flux in both directions is approximately zero,the reaction
is considered blocked. The function returns the IDs of the blocked reactions.

# INPUTS

- `myModel`:            A StandardModel that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

- `Tolerance`:          A small number that represents the level of error tolerance.
- `printLevel`:         Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `blocked_index`:      IDs of of blocked reactions.

# EXAMPLES

- Full input/output example
```julia
julia> blocked_index = find_blocked_reactions(myModel)
```

See also: `dataOfModel()`, `reversibility()`

"""

function find_blocked_reactions(myModel::StandardModel, Tolerance::Float64=1e-6, printLevel::Int=1)

    ## Export data from model

    S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel)

    ## Determine the reversibility of a reaction

    # Create an array of reaction IDs:
    Reaction_Ids = collect(1:n)
    irreversible_reactions_id, reversible_reactions_id = reversibility(lb, Reaction_Ids)
    n_irr = length(irreversible_reactions_id)
    n_rev = length(reversible_reactions_id)

    ## Homogenize the upper_bound and lower_bound of reactions

    printstyled("Homogenization:\n"; color=:cyan)

    # Set the maximum value for M:
    M = getM()

    ## Loop through each value in the array "lb" and "ub"

    for i in 1:n
        # If the lower bound is greater than zero, set it to zero:
        if lb[i] > 0
            lb[i] = 0
        end
        # If the upper bound is greater than zero, set it to M:
        if ub[i] > 0
            ub[i] = M
        end
        # If the lower bound is less than zero, set it to -M:
        if lb[i] < 0
            lb[i] = -M
        end
        # If the upper bound is less than zero, set it to zero:
        if ub[i] < 0
            ub[i] = 0
        end
    end

    ## Find Irreversible blocked reactions

    # Create empty arrays to hold the IDs of blocked and unblocked irreversible reactions:
    irreversible_blocked_reactions_id = []
    irreversible_unblocked_reactions_id = []

    # Create a new optimization model using the GLPK optimizer:
    model_irr = Model(GLPK.Optimizer)

    # Define the variable V for each reaction, with its lower and upper bounds:
    @variable(model_irr, lb[i] <= V[i = 1:n] <= ub[i])

    # Add a constraint to enforce that S*V = 0:
    @constraint(model_irr, S * V .== 0)

    ## Loop through each irreversible reaction in the irreversible reactions array

    for j in irreversible_reactions_id
        # Set the objective to maximize the flux of the current reaction:
        @objective(model_irr, Max, V[j])

        # Add a constraint to enforce that the flux of the current reaction is less than or equal to 1:
        @constraint(model_irr, c,  V[j] <= 1)

        # Optimize the model:
        optimize!(model_irr)

        # If the objective value is approximately zero, add the current reaction ID to the blocked reactions array:
        if isapprox(objective_value(model_irr), 0, atol = Tolerance)
            append!(irreversible_blocked_reactions_id, j)
        # Otherwise, add the current reaction ID to the unblocked reactions array:
        else
            append!(irreversible_unblocked_reactions_id, j)
        end

        # Remove the current constraint from the model:
        delete(model_irr, c)

        # Unregister the current constraint from the model:
        unregister(model_irr, :c)
    end

    ## Find Reversible blocked reactions

    # Create empty arrays to hold the IDs of blocked and unblocked reversible reactions:
    reversible_unblocked_reactions_id = []
    reversible_blocked_reactions_id = []

    # Create a new optimization model using the GLPK optimizer:
    model_rev = Model(GLPK.Optimizer)

    # Define the variable V for each reaction, with its lower and upper bounds:
    @variable(model_rev, lb[i] <= V[i = 1:n] <= ub[i])

    # Add a constraint to enforce that S*V = 0:
    @constraint(model_rev, S * V .== 0)

    # Loop through each reversible reaction in the reversible reactions array:
    for j in reversible_reactions_id

        ## Forward direction

        @objective(model_rev, Max, V[j]) # Set the objective to maximize the flux of the current reaction
        @constraint(model_rev, c1, V[j] <= 1) # Add a constraint to enforce that the flux of the current reaction is less than or equal to 1
        optimize!(model_rev) # Optimize the model
        opt_fwd = objective_value(model_rev) # Save the objective value
        t_fwd = termination_status(model_rev) # Save the termination status
        delete(model_rev, c1) # Remove the current constraint from the model
        unregister(model_rev, :c1) # Unregister the current constraint from the model

        ## Backward direction

        @objective(model_rev, Min, V[j]) # Set the objective to minimize the flux of the current reaction
        @constraint(model_rev, c2, V[j] >= -1) # Add a constraint to enforce that the flux of the current reaction is greater than or equal to -1
        optimize!(model_rev) # Optimize the model
        opt_back = objective_value(model_rev) # Save the objective value
        t_back = termination_status(model_rev) # Save the termination status
        delete(model_rev, c2) # Remove the current constraint from the model
        unregister(model_rev, :c2) # Unregister the current constraint from the model

        # If the objective values in both directions are approximately zero, add the current reaction ID to the blocked reactions array:
        if isapprox(opt_fwd, 0, atol = Tolerance) && isapprox(opt_back, 0, atol = Tolerance)
            append!(reversible_blocked_reactions_id, j)
        # Otherwise, add the current reaction ID to the unblocked reactions array:
        else
            append!(reversible_unblocked_reactions_id, j)
        end
    end

    ## Union the IDs of blocked irreversible and reversible reactions into a single array and Sort the array of blocked reaction IDs

    blocked_index = union(reversible_blocked_reactions_id, irreversible_blocked_reactions_id)
    blocked_index = sort(blocked_index)

    ## Print out results if requested

    if printLevel > 0
        printstyled("Consistency_Checking(TheNaiveApproch):\n"; color=:cyan)
        println("Number of Proccess : $(nprocs())")
        println("Number of Workers  : $(nworkers())")
        printstyled("Tolerance = $Tolerance\n"; color=:magenta)
        println("Number of irreversible blocked reactions : $(length(irreversible_blocked_reactions_id))")
        println("Number of reversible   blocked reactions : $(length(reversible_blocked_reactions_id))")
        println("Number of blocked reactions              : $(length(blocked_index))")
    end

    return blocked_index

end

end

#-----------------------------------------------------------------------------------------------------------------------------------------------------
