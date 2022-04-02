#-------------------------------------------------------------------------------------------

#=
    Purpose:    Preprocessing functions of metabolic networks
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
    setM(x)

Function that assigns a large value to M representing the concept of infinite boundary.

# INPUTS

- `x`:              a large number.

# OPTIONAL INPUTS

-

# OUTPUTS

-

# EXAMPLES

- Full input/output example
```julia
julia> setM(+1000.0)
```

"""

function setM(x)
    global M = x
    return
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
    check_duplicate_reaction(Reactions)

Function that examines metabolic networks to see if there is a repetitive reaction.

# INPUTS

- `Reactions`:      List of metabolic network reactions.

# OPTIONAL INPUTS

-

# OUTPUTS

- `true`:           There is a repetitive reaction.
- `false`:          There is no repetitive reaction.

# EXAMPLES

- Full input/output example
```julia
julia> check_duplicate = check_duplicate_reaction(Reactions)
```

"""

function check_duplicate_reaction(Reactions)
    n = length(Reactions)
    unique_reactions = unique!(Reactions)
    n_unique = length(unique_reactions)
    if n == n_unique
        return false
    else
        return true
    end
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
    reversibility_checking(lb, reversible_reactions_id)

Function that detects reversible reactions that are blocked in only one direction.

# INPUTS

- `lb`:                                   LowerBound Of Reactions.
- `reversible_reactions_id`:              Reversible reaction IDs.

# OPTIONAL INPUTS

-

# OUTPUTS

- `rev_blocked_fwd`:             Reversible reactions that are blocked in forward direction.
- `rev_blocked_back`:            Reversible reactions that are blocked in backward direction.

# EXAMPLES

- Full input/output example
```julia
julia> rev_blocked_fwd, rev_blocked_back = reversibility_checking(lb, reversible_reactions_id)
```

"""

function reversibility_checking(lb, reversible_reactions_id)

    n = length(lb)
    model = Model(GLPK.Optimizer)
    @variable(model, lb[i] <= V[i = 1:n] <= ub[i])
    @constraint(model, S * V .== 0)
    rev_blocked_fwd = []
    rev_blocked_back = []

    for j in reversible_reactions_id

    # The forward direction:

        @objective(model, Max, V[j])
        @constraint(model, V[j] <= 1)
        optimize!(model)
        opt_fwd = objective_value(model)

        if opt_fwd ≈ 0
             append!(rev_blocked_fwd, j)
        end

    # The backward direction:

        @objective(model, Min, V[j])
        @constraint(model, V[j] >= -1)
        optimize!(model)
        opt_back = objective_value(model)

        if opt_back ≈ 0
             append!(rev_blocked_back, j)
        end
     end
    return rev_blocked_fwd, rev_blocked_back
end

#-------------------------------------------------------------------------------------------

"""
    reversibility_correction(S, lb, ub, irreversible_reactions_id, reversible_reactions_id, rev_blocked_fwd, rev_blocked_back)

Function that modifies 3 sets:

    1) Remove rev_blocked_fwd and rev_blocked_back from reversible reactions list.
    2) Add rev_blocked_fwd and rev_blocked_back from irreversible reactions list.
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

- `S`:                                            Stoichiometric matrix.
- `lb`:                                           LowerBound Of Reactions.
- `ub`:                                           UpperBound of Reactions.
- `irreversible_reactions_id`:                    Irreversible reaction IDs.
- `corrected_reversible_reactions_id`:            Irreversible reaction IDs.

# EXAMPLES

- Full input/output example
```julia
julia> S, lb, ub, irreversible_reactions_id, corrected_reversible_reactions_id = reversibility_correction(S, lb, ub, irreversible_reactions_id, reversible_reactions_id, rev_blocked_fwd, rev_blocked_back)
```

"""

function reversibility_correction(S, lb, ub, irreversible_reactions_id, reversible_reactions_id, rev_blocked_fwd, rev_blocked_back)

    corrected_reversible_reactions_id = []

    # Forward

    for i in rev_blocked_fwd
    # Modify lb, ub:
        ub[i] = lb[i] * -1
        lb[i] = 0.0
    # Add to irreversible reactions list:
        append!(irreversible_reactions_id, i)
    # Modify S :
        S[:, i] .= S[:, i] * -1
    end

    # Backward

    for i in rev_blocked_back
    # Modify lb, ub:
        lb[i] = 0.0
    # Add to irreversible reactions list:
        append!(irreversible_reactions_id, i)
    end

    # Remove rev_blocked_fwd and rev_blocked_back from reversible reactions list

    set_reversible_reactions_id = Set(reversible_reactions_id)
    set_rev_blocked_fwd = Set(rev_blocked_fwd)
    set_rev_blocked_back = Set(rev_blocked_back)
    set_rev_blocked_onedirection = union(set_rev_blocked_fwd, set_rev_blocked_back)
    set_reversible_reactions_id = setdiff(set_reversible_reactions_id, set_rev_blocked_onedirection)

    for i in set_reversible_reactions_id
        append!(corrected_reversible_reactions_id, i)
    end
    return S, lb, ub, irreversible_reactions_id, corrected_reversible_reactions_id
end

export dataOfModel, setM, reversibility, check_duplicate_reaction, homogenization, reversibility_checking, reversibility_correction
