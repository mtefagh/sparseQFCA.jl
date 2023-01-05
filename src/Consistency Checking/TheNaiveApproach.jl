#-------------------------------------------------------------------------------------------

#=
    Purpose:    Finding blocked reactions in metabolic network
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology - Iran
    Date:       April 2022
=#

#-------------------------------------------------------------------------------------------

module TheNaiveApproach
export find_blocked_reactions

using GLPK, JuMP, COBREXA

include("../pre_processing.jl")

using .pre_processing

"""
    find_blocked_reactions(myModel)

A function that finds blocked reactions in metabolic network.

# INPUTS

- `myModel`:        A StandardModel that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

-

# OUTPUTS

- `blocked_names`:  Blocked reaction Names.

# EXAMPLES

- Full input/output example
```julia
julia> blocked_reactions = find_blocked_reactions(myModel)
```

See also: `dataOfModel()`, `reversibility()`, 'getTolerance()'

"""

function find_blocked_reactions(myModel::StandardModel)

    # Exporting data from model:

    S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel)

    # assigning a small value to Tolerance representing the level of error tolerance:

    Tolerance = getTolerance()

    # Determining the reversibility of a reaction:

    irreversible_reactions_id, reversible_reactions_id = reversibility(lb)

    # Homogenizing the upper_bound and lower_bound of reactions:

    M = 1000.0

    for i in 1:n
        if lb[i] > 0
            lb[i] = 0
        end
        if ub[i] > 0
            ub[i] = M
        end
        if lb[i] < 0
            lb[i] = -M
        end
        if ub[i] < 0
            ub[i] = 0
        end
    end

    # Irreversible blocked:

    irreversible_blocked_reactions_id = []
    irreversible_unblocked_reactions_id = []

    model_irr = Model(GLPK.Optimizer)
    @variable(model_irr, lb[i] <= V[i = 1:n] <= ub[i])
    @constraint(model_irr, S * V .== 0)
    for j in irreversible_reactions_id
        @objective(model_irr, Max, V[j])
        @constraint(model_irr, c,  V[j] <= 1)
        # @constraint(model, V[j] >= 0)
        optimize!(model_irr)
        if isapprox(objective_value(model_irr), 0, atol = Tolerance)
            append!(irreversible_blocked_reactions_id, j)
        else
            append!(irreversible_unblocked_reactions_id, j)
        end
        delete(model_irr, c)
        unregister(model_irr, :c)
    end

    # Reversible blocked:

    reversible_unblocked_reactions_id = []
    reversible_blocked_reactions_id = []

    model_rev = Model(GLPK.Optimizer)
    @variable(model_rev, lb[i] <= V[i = 1:n] <= ub[i])
    @constraint(model_rev, S * V .== 0)
    for j in reversible_reactions_id

        #Forward

        @objective(model_rev, Max, V[j])
        @constraint(model_rev, c1, V[j] <= 1)
        # @constraint(model, V[j] >= 0)
        optimize!(model_rev)
        opt_fwd = objective_value(model_rev)
        t_fwd = termination_status(model_rev)
        delete(model_rev, c1)
        unregister(model_rev, :c1)

        #Backward

        @objective(model_rev, Min, V[j])
        @constraint(model_rev, c2, V[j] >= -1)
        # @constraint(model, V[j] >= 0)
        optimize!(model_rev)
        opt_back = objective_value(model_rev)
        t_back = termination_status(model_rev)
        delete(model_rev, c2)
        unregister(model_rev, :c2)
        if isapprox(opt_fwd, 0, atol = Tolerance) && isapprox(opt_back, 0, atol = Tolerance)
            append!(reversible_blocked_reactions_id, j)
        else
            append!(reversible_unblocked_reactions_id, j)
        end
    end

    blocked_index = union(reversible_blocked_reactions_id, irreversible_blocked_reactions_id)
    blocked_index = sort(blocked_index)

    return blocked_index

end

end
