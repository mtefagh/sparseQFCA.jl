#-----------------------------------------------------------------------------------------------------------------------------------------------------
#=
    Purpose:    Metabolic Network Reductions based on Quantitative Flux Coupling Analysis and the concept of lossy compression in Information Theory.
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology
    Date:       May 2023
=#
#-----------------------------------------------------------------------------------------------------------------------------------------------------

module QuantomeReducer

export quantomeReducer

using COBREXA, SparseArrays, GLPK, JuMP, LinearAlgebra, Distributed, SharedArrays, Clarabel

import CDDLib

import JSONFBCModels: JSONFBCModel

import AbstractFBCModels as A

include("../Pre_Processing/Pre_processing.jl")
using .Pre_processing

include("../ConsistencyChecking/SwiftCC.jl")
using .SwiftCC

include("../QFCA/distributedQFCA.jl")
using .DistributedQFCA

"""
    quantomeReducer(model)

The function is designed to perform metabolic network reduction by removing blocked reactions, merge all the fully coupled reactions,
remove the eligible reactions by the DCE-induced reductions.It extracts relevant data, separates reversible and irreversible reactions,
corrects reversibility, and removes zero rows from the stoichiometric matrix. It processes flux coupling analysis to identify reaction
clusters and reactions to be removed. The function constructs a reduced metabolic network matrix and performs distributed optimization
for DCE-induced reductions. Finally, it generates information about the reduction process and returns the reduced metabolic network matrix.

# INPUTS

- `model`:                     A model that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

- `removing`:                  A boolean variable that indicates whether reactions should be removed from the network in the stages of determining coupling or not.
- `Tolerance`:                 A small number that represents the level of error tolerance.
- `printLevel`:                Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `reducedModel`               A reduced metabolic network.
- `A`                          A `n` x `ñ` matrix representing the coefficients betwee reactions of original networks and reactions of reduced network.
- `reduct-map`                 A dictionary to save the representatives of eliminations.

# EXAMPLES

- Full input/output example
```julia
julia> reducedModel, A, reduct_map = quantomeReducer(model)
```

See also: `dataOfModel()`, , `reversibility()`, `homogenization()`, `distributedReversibility_Correction()`, `distributedQFCA()`

"""

function quantomeReducer(model, removing::Bool=false, Tolerance::Float64=1e-6, OctuplePrecision::Bool=false, printLevel::Int=1)

    ## Extracte relevant data from input model

    S, Metabolites, Reactions, Genes, m, n, n_genes, lb, ub, c_vector = dataOfModel(model, printLevel)

    ## Find Biomass

    index_c = findfirst(x -> x == 1.0, c_vector)
    Biomass = Reactions[index_c]

    row_S, col_S = size(S)

    ## Ensure that the bounds of all reactions are homogenous

    lb, ub = homogenization(lb, ub, 0)

    ## Separate reactions into reversible and irreversible sets

    # Create an array of reaction IDs:
    Reaction_Ids = collect(1:n)
    irreversible_reactions_id, reversible_reactions_id = reversibility(lb, Reaction_Ids, printLevel)

    ## Create a new instance of the input model with homogenous bounds

    ModelObject_CC = Model_CC(S, Metabolites, Reactions, Genes, m, n, lb, ub)
    blocked_index, ν  = swiftCC(ModelObject_CC, Tolerance, OctuplePrecision, printLevel)
    blocked_index_rev = blocked_index ∩ reversible_reactions_id

    # Convert to Vector{Int64}:
    blocked_index_rev = convert(Vector{Int64}, blocked_index_rev)

    ## Correct Reversibility

    Dict_bounds = Dict()

    for i = 1:n
        Dict_bounds[Reactions[i]] = lb[i], ub[i]
    end

    ModelObject_Crrection = Model_Correction(S, Metabolites, Reactions, Genes, m, n, lb, ub, irreversible_reactions_id, reversible_reactions_id)

    S, lb, ub, irreversible_reactions_id, reversible_reactions_id = distributedReversibility_Correction(ModelObject_Crrection, blocked_index_rev, OctuplePrecision, printLevel)
    for i = 1:n
        Dict_bounds[Reactions[i]] = lb[i], ub[i]
    end

    row, col = size(S)
    model_Correction_Constructor(ModelObject_Crrection , S, Metabolites, Reactions, Genes, row, col, lb, ub, irreversible_reactions_id, reversible_reactions_id)

    ## Obtain blocked_index, fctable, Fc_Coefficients, and Dc_Coefficients

    row, col = size(S)
    ModelObject_QFCA = Model_QFCA(S, Metabolites, Reactions, Genes, row, col, lb, ub, irreversible_reactions_id, reversible_reactions_id)

    # Convert to Vector{Int64}:
    blocked_index = convert(Vector{Int64}, blocked_index)
    fctable, Fc_Coefficients, Dc_Coefficients = distributedQFCA(ModelObject_QFCA, blocked_index)

    # Get the dimensions of fctable:
    row, col = size(fctable)

    ## Count the number of reactions in each set

    n_irr = length(irreversible_reactions_id)
    n_rev = length(reversible_reactions_id)

    # Remove blocked reactions from Reaction_Ids:
    Reaction_Ids_noBlocked = setdiff(Reaction_Ids, blocked_index)

    ## Remove any reactions that cannot carry flux and is blocked

    Reactions_noBlocked = Reactions[setdiff(range(1, n), blocked_index)]
    lb_noBlocked = lb[setdiff(1:end, blocked_index)]
    ub_noBlocked = ub[setdiff(1:end, blocked_index)]
    S_noBlocked  = S[:, setdiff(1:end, blocked_index)]
    irreversible_reactions_id = setdiff(irreversible_reactions_id, blocked_index)
    reversible_reactions_id   = setdiff(reversible_reactions_id, blocked_index)

    ## Remove all rows from a given sparse matrix S that contain only zeros and the corresponding metabolites from the Metabolites array

    # Sort irreversible_reactions_id and reversible_reactions_id:
    irreversible_reactions_id = sort(irreversible_reactions_id)
    reversible_reactions_id = sort(reversible_reactions_id)

    # Initialize arrays:
    A_rows_original = Array{Int64}([])
    A_cols_reduced = Array{Int64}([])

    # Make copies of Reaction_Ids for later use:
    A_rows_original = copy(Reaction_Ids)
    A_cols_reduced = copy(Reaction_Ids)

    ## FC

    # Initialize dictionaries and lists:
    FC = Dict()
    FC_Final = Dict()
    FC_Coef = Dict()
    remove_list_FC = []
    c = 1

    ## Iterate over fctable to identify and store FC coefficients

    # Iterating over a range starting from 1 and ending at col:
    for i in range(1, col)
        # Nested loop iterating over a range starting from i+1 and ending at col:
        for j in range(i+1, col)
            # Checking conditions for equality and i not equal to j:
            if (fctable[i,j] == fctable[j,i] == 1.0) && (i != j)
                # Assigning tuple to FC_Coef[c]:
                FC_Coef[c] = Reaction_Ids_noBlocked[i], Reaction_Ids_noBlocked[j], Fc_Coefficients[i,j]
                # Assigning tuple to FC[c]:
                FC[c] = Reaction_Ids_noBlocked[i], Reaction_Ids_noBlocked[j]
                # Incrementing the counter c by 1:
                c = c + 1
            end
        end
    end

    ## Process FC dictionary to remove duplicates and create FC_Final dictionary

    s = 1
    # Iterating over keys of FC after sorting them:
    for key in sort(collect(keys(FC)))
        # Checking conditions for key comparison:
        if (key >= 2) && (FC[key][1] == FC[key-1][1])
            # Initializing an empty list called temp_list:
            temp_list = []
            # Initializing an empty list called delete_list:
            delete_list = []
            # Looping from 1 to key-1:
            for i = 1 : key-1
                # Checking condition for equality of FC elements:
                if FC[key][1] .== FC[i][1]
                    # Appending FC[i][2] to temp_list:
                    append!(temp_list, FC[i][2])
                    # Appending i to delete_list:
                    append!(delete_list, i)
                end
            end
            # Appending FC[key][2] to temp_list:
            append!(temp_list, FC[key][2])
            # Appending temp_list to remove_list_FC:
            append!(remove_list_FC, temp_list)
            # Assigning tuple to FC_Final[s]:
            FC_Final[s] = FC[key][1], temp_list
            # Incrementing the counter s by 1:
            s = s + 1
            # Looping over delete_list:
            for i in delete_list
                # Deleting element i from FC_Final:
                delete!(FC_Final, i)
            end
        else
            FC_Final[s] = FC[s]  # Assigning FC[s] to FC_Final[s]
            s = s + 1  # Incrementing the counter s by 1
        end
    end

    # Remove duplicate FC cluster members:
    remove_list_FC = unique(remove_list_FC)

    # Sort FC cluster members:
    remove_list_FC = sort(remove_list_FC)

    # Remove FC clusters that contain reactions in remove_list_FC from FC_Final dictionary:
    for key in sort(collect(keys(FC_Final)))
        if FC_Final[key][1] in remove_list_FC
            delete!(FC_Final, key)
        end
    end

    # Initialize arrays for FC representatives and cluster members:
    FC_representatives = Array{Int64}([])
    FC_cluster_members = Array{Int64}([])

    # Iterate through the items in the dictionary
    for (key, value) in FC_Final
        # Check if the second element of the tuple is an integer
        if isa(value[2], Integer)
            # Convert the integer to a Vector{Any}
            FC_Final[key] = (value[1], Any[value[2]])
        end
    end

    # Populate FC_representatives and FC_cluster_members arrays:
    for key in sort(collect(keys(FC_Final)))
        append!(FC_cluster_members, FC_Final[key][2])
        append!(FC_representatives, FC_Final[key][1])
        if index_c in FC_Final[key][2]
            index_c = FC_Final[key][1]
            Biomass = Reactions[index_c]
        end
    end

    # Create FC_Clusters dictionary:
    FC_Clusters = Dict()
    c = 1
    for key in sort(collect(keys(FC_Final)))
        FC_Clusters[c] = FC_Final[key]
        c += 1
    end

    ## DC

    # Initialize an empty array to store the IDs of reactions to be removed:
    remove_list_DC = Array{Int64}([])

    # Iterate over the rows:
    for i in range(1, row)
        # Iterate over the columns:
        for j in range(1, col)
            # Check if the value at position (i, j) in fctable is equal to 4.0:
            if 4.0 ∈ fctable[i, :]
                # If the condition is true, append the corresponding Reaction ID to remove_list_DC:
                append!(remove_list_DC, Reaction_Ids_noBlocked[i])
                # Exit the inner loop as the reaction has been found and added to remove_list_DC:
                break
            end
        end
    end

    # Sort the 'blocked_index', 'FC_cluster_members', and 'remove_list_DC' arrays:
    blocked_index = sort(blocked_index)
    FC_cluster_members = sort(FC_cluster_members)
    remove_list_DC = sort(remove_list_DC)
    Reaction_Ids_noBlocked = sort(Reaction_Ids_noBlocked)
    irreversible_reactions_id = sort(irreversible_reactions_id)
    reversible_reactions_id = sort(reversible_reactions_id)

    # Create the 'Eliminations' array by taking the union of 'blocked_index', 'remove_list_DC', and 'FC_cluster_members':

    Eliminations = union(blocked_index, remove_list_DC, FC_cluster_members)

    # Sort Eliminations:
    Eliminations = sort(Eliminations)

    # Update 'A_cols_reduced' by removing elements in 'Eliminations' from the range 1 to 'n' in 'A_cols_reduced':
    A_cols_reduced = A_cols_reduced[setdiff(range(1, n), Eliminations)]

    ## Matrix A

    # Create a shared array 'A' of size (A_rows_original, A_cols_reduced) with initial values set to false:
    n = length(A_rows_original)
    ñ = length(A_cols_reduced)
    A = SharedArray{Float64, 2}((n,ñ), init = false)

    ## Blocked

    # Iterate over indices from 1 to n:
    for i = 1:n
        if i in blocked_index
            # Set the entire row to 0.0 in 'A' for the indices present in 'blocked_index':
            A[i, :] .= 0.0
        end
    end

    ## I

    # Iterate over indices from 1 to ñ:
    for i = 1:ñ
        # Set the corresponding element to 1.0 in 'A' for each index in 'A_cols_reduced':
        A[A_cols_reduced[i], i] = 1.0
    end

    ## FC

    # Iterate over keys in 'FC_Coef' sorted in ascending order:
    for key in sort(collect(keys(FC_Coef)))
        if FC_Coef[key][1] in A_cols_reduced
            # Find the index in 'A_cols_reduced' where the first element of 'FC_Coef[key]' is present:
            index = findfirst(x -> x == FC_Coef[key][1], A_cols_reduced)
            # Set the corresponding element in 'A' to the third element of 'FC_Coef[key]':
            A[FC_Coef[key][2], index] = FC_Coef[key][3]
        end
    end

    ## DC

    DCE = Dict()
    counter = 1

    # Create a new optimization model:
    model_local = GenericModel{BigFloat}(Clarabel.Optimizer{BigFloat})
    settings = Clarabel.Settings(verbose = false, time_limit = 5)

    # Define variables λ, ν, and t:
    @variable(model_local, λ[1:n])
    @variable(model_local, ν[1:m])
    @variable(model_local, t)

    # Set the objective function to minimize t:
    @objective(model_local, Min, t)

    @constraint(model_local, [j in reversible_reactions_id], λ[j] == 0.0)

    # Perform distributed optimization for each reaction in remove_list_DC:
    @sync @distributed for i in remove_list_DC

        # save i as row number:
        row = i

        # Compute the linear combination for DCE:
        DCE_LinearCombination = setdiff(1:n, i)

        # Define constraints for the optimization problem:
        con1 = @constraint(model_local, [t; λ] in MOI.NormOneCone(1 + length(λ)))
        con2 = @constraint(model_local, λ == S' * ν)
        con3 = @constraint(model_local, [j in Eliminations ∩ DCE_LinearCombination], λ[j] == 0.0)
        con4 = @constraint(model_local, [j in DCE_LinearCombination], λ[j] >= 0.0)
        @constraint(model_local, con5, λ[i] == -1.0)

        # Solve the optimization problem:
        optimize!(model_local)

        # Incrementing the counter by 1:
        counter += 1

        # Get the values of λ:
        λ_vec = value.(λ)

        # Create an empty array for storing column indices:
        cols = Array{Int64}([])

        for i = 1:length(λ_vec)
            if λ_vec[i] > Tolerance
                index = findfirst(x -> x == i, A_cols_reduced)
                A[row, index] = λ_vec[i]
                if index_c == row
                    index_c = i
                    Biomass = Reactions[index_c]
                end
            end
        end

        ## Condition 1

        delete(model_local, con1) # Remove the current constraint from the model
        unregister(model_local, :con1) # Unregister the current constraint from the model

        ## Condition 2

        delete(model_local, con2) # Remove the current constraint from the model
        unregister(model_local, :con2) # Unregister the current constraint from the model

        ## Condition 3

        constraint_refs_con3 = [con3[i] for i in eachindex(con3)]
        for i in constraint_refs_con3
            delete(model_local, i)
            unregister(model_local, :i)
        end

        ## Condition 4

        constraint_refs_con4 = [con4[i] for i in eachindex(con4)]
        for i in constraint_refs_con4
            delete(model_local, i)
            unregister(model_local, :i)
        end

        ## Condition 5

        delete(model_local, con5) # Remove the current constraint from the model
        unregister(model_local, :con5) # Unregister the current constraint from the model

    end

    A = convert(Matrix{Float64}, A)
    row_A, col_A = size(A)

    ## Matrix S̃

    row_S, col_S = size(S)
    S̃ = S * A
    S̃, Metabolites_reduced, Metabolites_elimination = remove_zeroRows(S̃, Metabolites)
    S̃_row, S̃_col = size(S̃)

    # Copying A_cols_reduced and assigning it to R̃:
    Reactions_elimination = Reactions[setdiff(range(1, n), A_cols_reduced)]
    R̃ = Reactions[A_cols_reduced]

    # Copying Metabolites_reduced and assigning it to M̃:
    M̃ = copy(Metabolites_reduced)

    ## Reduction Map

    # Creating an empty dictionary called reduction_map:
    reduction_map = Dict()

    # Initializing the counter c to 1:
    c = 1

    # Iterating over each column of matrix A:
    for col in eachcol(A)
        # Converting the current column to a sparse vector:
        col = sparsevec(col)
        # Finding the non-zero indices and values in the sparse vector:
        nonzero_indices, nonzero_values = findnz(col)
        index = findfirst(x -> x == R̃[c], Reactions)
        # Assigning the tuple (R̃[c], nonzero_indices) to the key c in reduction_map:
        reduction_map[c] = index, nonzero_indices
        # Incrementing the counter c by 1:
        c += 1
    end

    ### Genes Modify

    for i in model.reactions
        for (key, value) in sort(reduction_map)
            if i.first == Reactions[value[1]]
                for j in value[2]
                    if !isnothing(i.second.gene_association_dnf) && !isnothing(model.reactions[Reactions[j]].gene_association_dnf)
                        i.second.gene_association_dnf = vcat(i.second.gene_association_dnf, model.reactions[Reactions[j]].gene_association_dnf)
                    end
                end
            end
        end
    end

    for i in model.reactions
        # Check if gene_association_dnf exists and is not nothing
        if !isnothing(i.second.gene_association_dnf)
            flattened = filter(!isempty, i.second.gene_association_dnf) |> Iterators.flatten

            seen = Set{String}()
            result = Vector{Vector{String}}()

            for gene in flattened
                if !(gene in seen)
                    push!(seen, gene)
                    push!(result, [gene])
                end
            end

            i.second.gene_association_dnf = result
        end

    end

    # lb & ub
    for i in model.reactions
        i.second.lower_bound = Dict_bounds[i.first][1]
        i.second.upper_bound = Dict_bounds[i.first][2]
    end

    # Reaction:

    filter!(pair -> !(pair.first in Reactions_elimination), model.reactions)

    filter!(pair -> !(pair.first in Metabolites_elimination), model.metabolites)

    for i in model.reactions
        for key in collect(keys(i.second.stoichiometry))
            delete!(i.second.stoichiometry, key)
        end
    end

    for i in model.reactions
        index_col = findfirst(x -> x == i.first, R̃)
        stoichiometry_vector = S̃[:,index_col]
        met = 1
        for c = 1:length(stoichiometry_vector)
            push!(i.second.stoichiometry, "$(M̃[met])" => stoichiometry_vector[c])
            met += 1
        end
    end


    for i in model.reactions
        if i.first == Biomass
            i.second.objective_coefficient = 1.0
        end
    end

    Genes_final = []
    for i in model.reactions
        if !isnothing(i.second.gene_association_dnf)
            append!(Genes_final, i.second.gene_association_dnf)
        end
    end

    Genes_final = unique(Genes_final)
    Genes_final = collect(Iterators.flatten(Genes_final))
    Genes_removal = setdiff(Genes, Genes_final)


    filter!(pair -> !(pair.first in Genes_removal), model.genes)

    ## Print out results if requested

    if printLevel > 0
        printstyled("Metabolic Network Reductions:\n"; color=:cyan)
        printstyled("Tolerance = $Tolerance\n"; color=:magenta)
        println("Original Network:")
        println("S           : $(row_S) x $(col_S)")
        println("Genes       : $(length(Genes))")
        println("Metabolites : $(m)")
        println("Reactions   : $(n)")
        println("Reduced Network:")
        println("S           : $(S̃_row) x $(S̃_col)")
        #println("Genes       : $(length(Genes_reduced))")
        println("Metabolites : $(length(Metabolites_reduced))")
        println("Reactions   : $(length(R̃))")
        println("A matrix    : $(row_A) x $(col_A)")
    end

    return model
end

end

#-----------------------------------------------------------------------------------------------------------------------------------------------------
