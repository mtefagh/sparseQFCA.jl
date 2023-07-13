#-------------------------------------------------------------------------------------------
#=
    Purpose:    Metabolic Network Reductions based on Quantitative Flux Coupling Analysis and the concept of lossy compression in Information Theory.
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology
    Date:       May 2023
=#
#-------------------------------------------------------------------------------------------

module Reduction
export reduction

include("../QFCA/distributedQFCA.jl")
include("../Data Processing/Pre_processing.jl")

using .Pre_processing, .DistributedQFCA, COBREXA, SparseArrays, GLPK, JuMP, LinearAlgebra, Distributed, SharedArrays, SparseArrays

"""
    reduction(myModel)

The function is designed to perform metabolic network reduction by removing blocked reactions, merge all the fully coupled reactions,
remove the eligible reactions by the DCE-induced reductions.It extracts relevant data, separates reversible and irreversible reactions,
corrects reversibility, and removes zero rows from the stoichiometric matrix. It processes flux coupling analysis to identify reaction clusters
and reactions to be removed. The function constructs a reduced metabolic network matrix and performs distributed optimization for DCE-induced reductions.
Finally, it generates information about the reduction process and returns the reduced metabolic network matrix.

# INPUTS

- `myModel`:                   A model that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

- `removing`:                  A boolean variable that indicates whether reactions should be removed from the network in the stages of determining coupling or not.
- `Tolerance`:                 A small number that represents the level of error tolerance.
- `printLevel`:                Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `A`                          A `n` x `ñ` matrix representing the reduced metabolic network matrix.

# EXAMPLES

- Full input/output example
```julia
julia> A = reduction(myModel)
```

See also: `dataOfModel()`, , `reversibility()`, `homogenization()`, `distributedReversibility_Correction()`, `distributedQFCA()`

"""

function reduction(myModel::StandardModel, removing::Bool=false, Tolerance::Float64=1e-6, printLevel::Int=1)

    ## Extracte relevant data from input model

    S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel, printLevel)

    ## Obtain blocked_index, fctable, Fc_Coefficients, and Dc_Coefficients

    blocked_index, fctable, Fc_Coefficients, Dc_Coefficients = distributedQFCA(myModel, removing, Tolerance, 0)

    # Get the dimensions of fctable:
    row, col = size(fctable)

    ## Ensure that the bounds of all reactions are homogenous

    lb, ub = homogenization(lb, ub, printLevel)

    ## Separate reactions into reversible and irreversible sets

    # Create an array of reaction IDs:
    Reaction_Ids = collect(1:n)
    irreversible_reactions_id, reversible_reactions_id = reversibility(lb, Reaction_Ids, printLevel)

    ## Correct Reversibility

    S, lb, ub, irreversible_reactions_id, reversible_reactions_id = distributedReversibility_Correction(S, lb, ub, irreversible_reactions_id, reversible_reactions_id, printLevel)
    row_S, col_S = size(S)

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

    S_noBlocked, Metabolites = remove_zeroRows(S_noBlocked,Metabolites)

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

    # Populate FC_representatives and FC_cluster_members arrays:
    for key in sort(collect(keys(FC_Final)))
        append!(FC_cluster_members, FC_Final[key][2])
        append!(FC_representatives, FC_Final[key][1])
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

    # Iterate over the rows
    for i in range(1, row)
        # Iterate over the columns
        for j in range(1, col)
            # Check if the value at position (i, j) in fctable is equal to 4.0:
            if (fctable[i, j] == 4.0 || fctable[i, j] == 2.0 )
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

    # Perform distributed optimization for each reaction in remove_list_DC:
    @sync @distributed for i in remove_list_DC

        # Compute the linear combination for DCE:
        DCE_LinearCombination = setdiff(1:n, i)

        # Create a new optimization model:
        model = Model(GLPK.Optimizer)

        # Define variables λ, dualVar, and t:
        @variable(model, λ[1:n])
        @variable(model, dualVar[1:m])
        @variable(model, t)

        # Set the objective function to minimize t:
        @objective(model, Min, t)

        # Define constraints for the optimization problem:
        @constraint(model, [t; λ] in MOI.NormOneCone(1 + length(λ)))
        @constraint(model, λ == S' * dualVar)
        @constraint(model, [j in Eliminations ∩ DCE_LinearCombination], λ[j] == 0.0)
        #@constraint(model, [j in DCE_LinearCombination], λ[j] >= 0.0)
        @constraint(model, λ[i] == -1.0)
        @constraint(model, [j in reversible_reactions_id], λ[j] == 0.0)

        # Solve the optimization problem:
        optimize!(model)

        # Incrementing the counter by 1:
        counter += 1

        # Get the values of λ:
        λ_vec = value.(λ)

        # Create an empty array for storing column indices:
        cols = Array{Int64}([])

        # Find the index of the negative value and collect the indices of positive values:
        for i = 1:length(λ_vec)
            if λ_vec[i] < 0.0
                row = i
            end
            if λ_vec[i] > 0.0
                append!(cols, i)
            end
        end

        # Update the corresponding elements of A with the values of λ:
        for i = 1:length(cols)
            A[row, i] = λ_vec[cols[i]]
        end
    end

    A = convert(Matrix{Float64}, A)
    row_A, col_A = size(A)

    ## Matrix S̃

    S̃ = S * A
    row_S̃, col_S̃ = size(S̃)

    ## R̃

    # Copying A_cols_reduced and assigning it to R̃:
    R̃ = copy(A_cols_reduced)

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
        # Assigning the tuple (R̃[c], nonzero_indices) to the key c in reduction_map:
        reduction_map[c] = R̃[c], nonzero_indices
        # Incrementing the counter c by 1:
        c += 1
    end

    ## Print out results if requested

    if printLevel > 0
        printstyled("Metabolic Network Reductions:\n"; color=:cyan)
        printstyled("Tolerance = $Tolerance\n"; color=:magenta)
        println("S : $(row_S) x $(col_S)")
        println("R : $(n)")
        println("S̃ : $(row_S̃) x $(col_S̃)")
        println("R̃ : $(ñ)")
        println("A : $(row_A) x $(col_A)")
    end

    return A
end

end
