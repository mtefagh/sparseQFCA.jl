#-------------------------------------------------------------------------------------------
#=
    Purpose:
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology
    Date:       May 2023
=#
#-------------------------------------------------------------------------------------------

module Reduction
export reduction

include("../QFCA/DistributedQFCA.jl")
include("../Data Processing/Pre_processing.jl")

using .Pre_processing, .DistributedQFCA, COBREXA, SparseArrays, GLPK, JuMP, LinearAlgebra, Distributed, SharedArrays, SparseArrays

"""
    reduction(myModel)

# INPUTS

- `myModel`:                   A model that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

- `removing`:                  A boolean variable that indicates whether reactions should be removed from the network in the stages of determining coupling or not.
- `Tolerance`:                 A small number that represents the level of error tolerance.
- `printLevel`:                Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `A`

# EXAMPLES

- Full input/output example
```julia
julia>
```

See also: `dataOfModel()`, , `reversibility()`, `homogenization()`, `distributedReversibility_Correction()`, `distributedQFCA()`

"""

function reduction(myModel::StandardModel, removing::Bool=false, Tolerance::Float64=1e-6, printLevel::Int=1)

    ## Extract data from the model

    S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel, 0)

    ## Homogenize lower and upper bounds

    lb, ub = homogenization(lb, ub, 0)

    ## Determine irreversible and reversible reactions

    irreversible_reactions_id, reversible_reactions_id = reversibility(lb, 0)

    ## Remove zero rows from S matrix

    S, Metabolites = remove_zeroRows(S, Metabolites)
    m = length(Metabolites)

    ## Perform distributed reversibility correction

    S, lb, ub, irreversible_reactions_id, reversible_reactions_id = distributedReversibility_Correction(S, lb, ub, irreversible_reactions_id, reversible_reactions_id, 0)

    ## Obtain blocked_index, fctable, Fc_Coefficients, and Dc_Coefficients

    blocked_index, fctable, Fc_Coefficients, Dc_Coefficients = distributedQFCA(myModel, removing, Tolerance, 0)

    # Get the dimensions of fctable:
    row, col = size(fctable)

    # Initialize arrays:
    A_rows_original = Array{Int64}([])
    A_cols_reduced = Array{Int64}([])

    # Create an array of reaction IDs:
    Reaction_Ids = collect(1:n)

    # Remove blocked reactions from Reaction_Ids:
    Reaction_Ids_noBlocked = setdiff(Reaction_Ids, blocked_index)

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

    remove_list_DC = Array{Int64}([])

    # Iterate over rows in the range from 1 to 'row':
    for i in range(1, row)
        DC_row = Dc_Coefficients[i, :]
        DC_row = sparsevec(DC_row)
        # Check if the number of non-zero elements in the row is greater than or equal to 2:
        if nnz(DC_row) >= 2
            for j = 1:length(DC_row)
                # Check if the element is approximately zero based on the given tolerance:
                if isapprox(DC_row[j], 0, atol=Tolerance)
                    Dc_Coefficients[i, j] = 0.0
                end
            end
            DC_row = Dc_Coefficients[i, :]
            DC_row = sparsevec(DC_row)
            append!(remove_list_DC, i)
        end
    end

    # Sort the 'blocked_index', 'FC_cluster_members', and 'remove_list_DC' arrays:
    blocked_index = sort(blocked_index)
    FC_cluster_members = sort(FC_cluster_members)
    remove_list_DC = sort(remove_list_DC)

    # Create the 'Eliminations' array by taking the union of 'blocked_index', 'remove_list_DC', and 'FC_cluster_members':
    Eliminations = union(blocked_index, remove_list_DC, FC_cluster_members)

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

    # Perform distributed optimization for each reaction in remove_list_DC:
    @sync @distributed for i in remove_list_DC

        # Compute the linear combination for DCE:
        DCE_LinearCombination = setdiff(1:n, i)

        # Create a new optimization model:
        model = Model(GLPK.Optimizer)

        # Define variables lambda, dualVar, and t:
        @variable(model, lambda[1:n])
        @variable(model, dualVar[1:m])
        @variable(model, t)

        # Set the objective function to minimize t:
        @objective(model, Min, t)

        # Define constraints for the optimization problem:
        c1 = @constraint(model, [t; lambda] in MOI.NormOneCone(1 + length(lambda)))
        c2 = @constraint(model, lambda == S' * dualVar)
        c3 = @constraint(model, [j in Eliminations ∩ DCE_LinearCombination], lambda[j] == 0.0)
        c4 = @constraint(model, [j in DCE_LinearCombination], lambda[j] >= 0.0)
        c5 = @constraint(model, lambda[i] == -1.0)
        c6 = @constraint(model, [j in reversible_reactions_id], lambda[j] == 0.0)

        # Solve the optimization problem:
        optimize!(model)

        # Get the values of lambda:
        lambda_vec = value.(lambda)

        # Create an empty array for storing column indices:
        cols = Array{Int64}([])

        # Find the index of the negative value and collect the indices of positive values:
        for i = 1:length(lambda_vec)
            if lambda_vec[i] < 0.0
                row = i
            end
            if lambda_vec[i] > 0.0
                append!(cols, i)
            end
        end

        # Update the corresponding elements of A with the values of lambda:
        for i = 1:length(cols)
            A[row, i] = lambda_vec[cols[i]]
        end
    end

    A = convert(Matrix{Float64}, A)

    ## S̃ Matrix

    S̃ = S * A

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

    return A
end

end
