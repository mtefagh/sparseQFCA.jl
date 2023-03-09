#-------------------------------------------------------------------------------------------
#=
    Purpose:    Preprocessing functions of metabolic networks analysis
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology - Iran
    Date:       April 2022
=#
#-------------------------------------------------------------------------------------------

module pre_processing
export MyModel, myModel_Constructor, dataOfModel, getM, getTolerance, reversibility, check_duplicate_reaction, homogenization, reversibility_checking, reversibility_correction

using Distributed

@everywhere using GLPK, JuMP, COBREXA, SparseArrays

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

@everywhere mutable struct MyModel
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

@everywhere function myModel_Constructor(ModelObject::MyModel, S::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}, Metabolites::Array{String,1}, Reactions::Array{String,1},
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
    dataOfModel(myModel)

A function that extracts essential data from StandardModel that has been built using COBREXA's `load_model` function.

# INPUTS

- `myModel`:        A StandardModel that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

-

# OUTPUTS

- `S`:              Stoichiometric matrix.
- `Metabolites`:    List of metabolic network metabolites.
- `Reactions`:      List of metabolic network reactions.
- `Genes`:          List of metabolic network reactions.
- `m`:              Number of rows of stoichiometric matrix.
- `n`:              Number of columns of stoichiometric matrix.
- `lb`:             Lower bound vector (n x 1)
- `ub`:             Upper bound vector (n x 1)

# EXAMPLES

- Full input/output example
```julia
julia> S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel)
```

"""

function dataOfModel(myModel::StandardModel)

    # Extracting Data:
    S = stoichiometry(myModel)
    Metabolites = metabolites(myModel)
    Reactions = reactions(myModel)
    Genes = genes(myModel)
    m = length(metabolites(myModel))
    n = length(reactions(myModel))
    lb = lower_bounds(myModel)
    ub = upper_bounds(myModel)

    # Sorting Reactions:
    p = sortperm(Reactions)
    Reactions = Reactions[p]
    lb = lb[p]
    ub = ub[p]
    S = S[:,p]
    return S, Metabolites, Reactions, Genes, m, n, lb, ub
end

#-------------------------------------------------------------------------------------------

"""
    getM()

A function that returns a large value to set M representing the concept of infinite boundary.

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

@everywhere function getM()
    f = open("../config/ConfigFile.txt", "r")
    c = 0
    for lines in readlines(f)
        c = c + 1
        if c == 1
            M = lines
            M = parse(Float64, M)
            return M
        else
            return
        end
    end
end

#-------------------------------------------------------------------------------------------

"""
    getTolerance()

A function that returns a small value to set Tolerance representing the level of error tolerance.

# INPUTS

-

# OPTIONAL INPUTS

- `Tolerance`:              a small number.

# OUTPUTS

-

# EXAMPLES

- Full input/output example
```julia
julia> Tolerance = getTolerance()
```

"""

@everywhere function getTolerance()
    f = open("../config/ConfigFile.txt", "r")
    c = 0
    for lines in readlines(f)
        c = c + 1
        if c == 1
            continue
        else
            Tolerance = lines
            Tolerance = parse(Float64, Tolerance)
            return Tolerance
        end
    end
end

#-------------------------------------------------------------------------------------------

"""
    reversibility(lb)

A function that determines the reversibility of a reaction from the lower_bound of a reaction.

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

See also: `dataOfModel()`

"""

@everywhere function reversibility(lb::Array{Float64,1})
    n = length(lb)
    irreversible_reactions_id = Array{Int64}([])
    reversible_reactions_id = Array{Int64}([])
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

A function that examines metabolic networks to see if there is a repetitive reaction.

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

See also: `dataOfModel()`

"""

function check_duplicate_reaction(Reactions::Array{String,1})
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

A function that homogenizes the upper_bound and lower_bound of reactions.

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

See also: `dataOfModel()`, 'getM()'

"""

@everywhere function homogenization(lb::Array{Float64,1}, ub::Array{Float64,1})
    n = length(lb)
    # Set a large number for M:
    M = getM()

    # If the lower bound is greater than zero, set it to zero
    lb[lb .> 0] .= 0
    # If the upper bound is greater than zero, set it to the constant M
    ub[ub .> 0] .= M
    # If the lower bound is less than zero, set it to the negative of the constant M
    lb[lb .< 0] .= -M
    # If the upper bound is less than zero, set it to zero
    ub[ub .< 0] .= 0

    return lb,ub
end

#-------------------------------------------------------------------------------------------

"""
    reversibility_checking(S, lb, ub, reversible_reactions_id)

A function that detects reversible reactions that are blocked in only one direction.

# INPUTS

- `S`:                                    Stoichiometric matrix.
- `lb`:                                   LowerBound Of Reactions.
- `ub`:                                   UpperBound of Reactions.
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

See also: `dataOfModel()`, `reversibility()`, 'getTolerance()'

"""

function reversibility_checking(S::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}, lb::Array{Float64,1}, ub::Array{Float64,1}, reversible_reactions_id::Vector{Int64})

    n = length(lb)
    Tolerance = getTolerance()
    model = Model(GLPK.Optimizer)
    @variable(model, lb[i] <= V[i = 1:n] <= ub[i])
    @constraint(model, S * V .== 0)
    rev_blocked_fwd = Array{Int64}([])
    rev_blocked_back = Array{Int64}([])

## Detection Loop ...

    for j in reversible_reactions_id

    # The Forward Direction:

        @objective(model, Max, V[j])
        @constraint(model, c1, V[j] <= 1)
        optimize!(model)
        opt_fwd = objective_value(model)

        if isapprox(opt_fwd, 0, atol=Tolerance)
             append!(rev_blocked_fwd, j)
        end
        delete(model, c1)
        unregister(model, :c1)

    # The Backward Direction:

        @objective(model, Min, V[j])
        @constraint(model, c2, V[j] >= -1)
        optimize!(model)
        opt_back = objective_value(model)

        if isapprox(opt_back, 0, atol=Tolerance)
             append!(rev_blocked_back, j)
        end
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

See also: `dataOfModel()`, `homogenization()`, `reversibility()`, reversibility_checking()

"""

function reversibility_correction(S::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}, lb::Array{Float64,1}, ub::Array{Float64,1}, irreversible_reactions_id::Vector{Int64},
                                  reversible_reactions_id::Vector{Int64}, rev_blocked_fwd::Vector{Int64}, rev_blocked_back::Vector{Int64})

    corrected_reversible_reactions_id = Array{Int64}([])

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

    # Remove rev_blocked_fwd and rev_blocked_back from reversible reactions list:

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

end
