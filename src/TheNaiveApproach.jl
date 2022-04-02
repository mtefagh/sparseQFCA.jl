#-------------------------------------------------------------------------------------------

#=
    Purpose:    Finding blocked reactions in metabolic network
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology - Iran
    Date:       April 2022
=#

#-------------------------------------------------------------------------------------------

"""
    dataOfModel(myModel)

Function that exports essential data from StandardModel that has been built using COBREXA's `load_model` function.

# INPUTS

- `myModel`:        A model that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

-

# OUTPUTS

- `S`:              Stoichiometric matrix.
- `Metabolites`:    List of metabolic network metabolites.
- `Reactions`:      List of metabolic network reactions.
- `Genes`:          List of metabolic network reactions.
- `m`:              Number of rows of stoichiometric matrix.
- `n`:              Number of columns of stoichiometric matrix.
- `lb`:             LowerBound Of Reactions.
- `ub`:             UpperBound of Reactions.

# EXAMPLES

- Full input/output example
```julia
julia> S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel)
```

"""

function dataOfModel(myModel)
    S = stoichiometry(myModel)
    Metabolites = metabolites(myModel)
    Reactions = reactions(myModel)
    Genes = genes(myModel)
    m = length(metabolites(myModel))
    n = length(reactions(myModel))
    lb = lower_bounds(myModel)
    ub = upper_bounds(myModel)
    return S, Metabolites, Reactions, Genes, m, n, lb, ub
end

#-------------------------------------------------------------------------------------------

"""
    reversibility(lb)

Function that determines the reversibility of a reaction from the lower_bound of a reaction.

# INPUTS

- `lb`:             LowerBound Of Reactions.

# OPTIONAL INPUTS

-

# OUTPUTS

- `irreversible_reactions_id`:            Irreversible reaction IDs.
- `reversible_reactions_id`:              Reversible reaction IDs.

# EXAMPLES

- Full input/output example
```julia
julia> irreversible_reactions_id, reversible_reactions_id = reversibility(lb)
```

"""

function reversibility(lb)
    n = length(lb)
    irreversible_reactions_id = []
    reversible_reactions_id = []
    for i in 1:n
        if lb[i] >= 0
            append!(irreversible_reactions_id, i)
        else
            append!(reversible_reactions_id, i)
        end
    end
    return irreversible_reactions_id, reversible_reactions_id
end

#-------------------------------------------------------------------------------------------

"""
    homogenization(lb,ub)

Function that homogenizes the upper_bound and lower_bound of reactions.

# INPUTS

- `lb`:             LowerBound Of Reactions.
- `ub`:             UpperBound of Reactions.

# OPTIONAL INPUTS

-

# OUTPUTS

- `lb`:             LowerBound Of Reactions has become homogenous.
- `ub`:             UpperBound of Reactions has become homogenous.

# EXAMPLES

- Full input/output example
```julia
julia> lb,ub = homogenization(lb,ub)
```

"""

function homogenization(lb,ub)
    n = length(lb)
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
    return lb,ub
end

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

"""

function find_blocked_reactions(myModel)

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
