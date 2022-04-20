#-------------------------------------------------------------------------------------------

#=
    Purpose:    Finding blocked reactions in metabolic network
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology - Iran
    Date:       April 2022
=#

#-------------------------------------------------------------------------------------------

"""
    find_blocked_reactions(myModel)

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
julia> blocked_reactions = find_blocked_reactions(myModel)
```

See also: `dataOfModel()`, `homogenization()`, `reversibility()`

"""

function find_blocked_reactions(myModel)
    
    # Exporting data from model
    
    S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel)
    
    # Determining the reversibility of a reaction 
    
    irreversible_reactions_id = []
    reversible_reactions_id = []
    irreversible_reactions_id, reversible_reactions_id = reversibility(lb)
    
    # Homogenizing the upper_bound and lower_bound of reactions
    
    lb,ub = homogenization(lb,ub)

    # Irreversible blocked

    irreversible_blocked_reactions_id = []
    irreversible_unblocked_reactions_id = []

    model_irr = Model(GLPK.Optimizer)
    @variable(model_irr, lb[i] <= V[i = 1:n] <= ub[i])
    @constraint(model_irr, S * V .== 0)
    for j in irreversible_reactions_id
        @objective(model_irr, Max, V[j])
        @constraint(model_irr,  V[j] <= 1)
        # @constraint(model, V[j] >= 0)
        optimize!(model_irr)
        if isapprox(objective_value(model_irr), 0, atol=1e-8)
            append!(irreversible_blocked_reactions_id, j)
        else
            append!(irreversible_unblocked_reactions_id, j)
        end
        # delete(model, c)
        # unregister(model, :c)
    end

    # Reversible blocked

    reversible_unblocked_reactions_id = []
    reversible_blocked_reactions_id = []

    model_rev = Model(GLPK.Optimizer)
    @variable(model_rev, lb[i] <= V[i = 1:n] <= ub[i])
    @constraint(model_rev, S * V .== 0)
    for j in reversible_reactions_id

        #Forward

        @objective(model_rev, Max, V[j])
        @constraint(model_rev, V[j] <= 1)
        # @constraint(model, V[j] >= 0)
        optimize!(model_rev)
        opt_fwd = objective_value(model_rev)
        t_fwd = termination_status(model_rev)
        # delete(model, c1)
        # unregister(model, :c1)

        #Backward

        @objective(model_rev, Min, V[j])
        @constraint(model_rev, V[j] >= -1)
        # @constraint(model, V[j] >= 0)
        optimize!(model_rev)
        opt_back = objective_value(model_rev)
        t_back = termination_status(model_rev)
        # delete(model, c2)
        # unregister(model, :c2)
        if isapprox(opt_fwd, 0, atol=1e-8) && isapprox(opt_back, 0, atol=1e-8)
            append!(reversible_blocked_reactions_id, j)
        else
            append!(reversible_unblocked_reactions_id, j)
        end
    end

    blocked_reactions = union(reversible_blocked_reactions_id, irreversible_blocked_reactions_id)
    blocked_names = []

    for i in blocked_reactions
        r_name = reactions(myModel)[i]
        push!(blocked_names, r_name)
    end

    blocked_names = sort(blocked_names)

    return blocked_names

end

export dataOfModel, reversibility, homogenization, find_blocked_reactions
