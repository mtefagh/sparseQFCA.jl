module pre_processing

export dataOfModel,setM,reversibility,check_duplicate_reaction,
                homogenization,reversibility_correctionremove_wrong_reversible

using COBREXA, JuMP, GLPK

#=
    dataOfModel
        input : model
        output : S, Metabolites, Reactions, Genes,
                #Row of S (m) = #Metabolites,
                #Culomn of S(n) = #Reactions,
                LowerBound Of Reactions(lb),
                UpperBounds of Reactions(ub)
=#

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

#=
    setM
        input : A large number
        output : M has been set
=#

function setM(x)
    global M = x
    return
end

#=
    Reversibility
        input : Reactions,n,lb
        output : irreversible_reactions_id, reversible_reactions_id
=#

function reversibility(Reactions,n,lb)
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

#=
    check_duplicate_reaction
        input : Reactions
        output :
        True : There are a number of repetitive reactions
        Flase : There is no repetitive Reaction
=#

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

#=
    homogenization
        input : lb,ub
        output : modified lb,ub
=#

function homogenization(lb,ub)
    irreversible_reactions_id = []
    reversible_reactions_id = []
    irreversible_reactions_id,reversible_reactions_id = reversibility(myModel)
    for i in irreversible_reactions_id
        if lb[i] > 0
            lb[i] = 0
        end
        if ub[i] > 0
            ub[i] = M
        end
    end
    for i in reversible_reactions_id
        if lb[i] < 0
            lb[i] = -M
        end
        if ub[i] > 0
            ub[i] = M
        end
    end
    return lb,ub
end

#=
    reversibility_correction
        input : S, irreversible_reactions_id, reversible_reactions_id, lb, ub
        output : Modified S, lb, ub
                 reversible_reactions_id : Unchanged
        irreversible_reactions_id : add Reversible reactions that are blocked in Forward or Backward Direction
                 rev_blocked_fwd : Reversible reactions that are blocked in Forward Direction
                 rev_blocked_back : Reversible reactions that are blocked in Backward Direction
=#

function reversibility_correction(S, irreversible_reactions_id, reversible_reactions_id, lb, ub)
    model = Model(GLPK.Optimizer)
    @variable(model, lb[i] <= V[i = 1:n] <= ub[i])
    @constraint(model, S * V .== 0)
    rev_blocked_fwd = []
    rev_blocked_back = []

## Detection Loop ...

    for j in reversible_reactions_id

    # The Forward Direction :

        @objective(model, Max, V[j])
        @constraint(model, V[j] <= 1)
    #   @constraint(model, V[j] >= 0)
        optimize!(model)
        opt_fwd = objective_value(model)

    # The Backward Direction :

        @objective(model, Min, V[j])
        @constraint(model, V[j] >= -1)
    #   @constraint(model, V[j] >= 0)
        optimize!(model)
        opt_back = objective_value(model)
        if opt_fwd ≈ 0
             append!(rev_blocked_fwd, j)
        end
        if opt_back ≈ 0
             append!(rev_blocked_back, j)
        end
     end

## Coorection Loop ...

    # Fowrard

     for i in rev_blocked_fwd
    # Update lb,ub :
        lb[i] = 0
        ub[i] = ub[i] * -1
    # Add to Irreversible Reactions :
        append!(irreversible_reactions_id, i)
    # Modify S :
        S[:, i] .= S[:, i] * -1
    # Remove From Reversible List :
        # index = splice!(reversible_reactions_id, i)
        # deleteat!(reversible_reactions_id,index)
      end

    # Backward

      for i in rev_blocked_back
    # Update lb,ub :
        lb[i] = 0
    # Add to Irreversible Reactions :
        append!(irreversible_reactions_id, i)
    # Remove From Reversible List :
         #  index = splice!(reversible_reactions_id, i)
         #  deleteat!(reversible_reactions_id,index)
      end
return S, irreversible_reactions_id, reversible_reactions_id, lb, ub, rev_blocked_fwd,rev_blocked_back
end

#=
    remove_wrong_reversible
        input  :  reversible_reactions_id,rev_blocked_fwd,rev_blocked_back
        output :  reversible_reactions_id : Corrected

=#

function remove_wrong_reversible(reversible_reactions_id,rev_blocked_fwd,rev_blocked_back)
    set_reversible_reactions_id = Set(reversible_reactions_id)
    set_rev_blocked_fwd = Set(rev_blocked_fwd)
    set_rev_blocked_back = Set(rev_blocked_back)
    union!(set_rev_blocked_fwd, set_rev_blocked_back)
    setdiff!(set_reversible_reactions_id, set_rev_blocked_fwd)
    reversible_reactions_id = list(set_reversible_reactions_id)
    return reversible_reactions_id
end

end
