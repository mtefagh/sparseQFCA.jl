#-----------------------------------------------------------------------------------------------------------------------------------------------------
#=
    Purpose:    Metabolic Network Reductions based on Quantitative Flux Coupling Analysis and the concept of lossy compression in Information Theory.
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology
    Date:       May 2023
=#
#-----------------------------------------------------------------------------------------------------------------------------------------------------

module QuantomeReducer

export quantomeReducer

using COBREXA, SparseArrays, JuMP, LinearAlgebra, Distributed, SharedArrays

import CDDLib
import Clarabel

import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel: Model
import AbstractFBCModels.CanonicalModel: Reaction, Metabolite, Gene, Coupling

include("../Pre_Processing/Solve.jl")
using .Solve

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

- `model`:                     A CanonicalModel that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

- `SolverName`:                Name of the solver(default: HiGHS).
- `OctuplePrecision`:          A flag(default: false) indicating whether octuple precision should be used when solving linear programs.
- `removing`:                  A flag controlling whether reactions should be filtered out during the coupling determination phase of network analysis.
- `Tolerance`:                 A small number that represents the level of error tolerance.
- `printLevel`:                Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `reducedModel`               A reduced metabolic network.

# EXAMPLES

- Full input/output example
```julia
julia> reducedModel = quantomeReducer(model)
```

See also: `dataOfModel()`, , `reversibility()`, `homogenization()`, `distributedReversibility_Correction()`, `distributedQFCA()`

"""

function quantomeReducer(model, SolverName::String="HiGHS", OctuplePrecision::Bool=false, removing::Bool=false, Tolerance::Float64=1e-6, printLevel::Int=1)

    ## Extracte relevant data from input model

    S, Metabolites, Reactions, Genes, m, n, n_genes, lb, ub, c_vector = dataOfModel(model, printLevel)

    # Find the index of the first occurrence where the element in c_vector is equal to 1.0:
    index_c = findfirst(x -> x == 1.0, c_vector)

    # Use the found index to retrieve the corresponding element from the Reactions array:
    Biomass = Reactions[index_c]

    row_S, col_S = size(S)

    ## Ensure that the bounds of all reactions are homogenous

    lb, ub = homogenization(lb, ub)

    ## Separate reactions into reversible and irreversible sets

    # Create an array of reaction IDs:
    Reaction_Ids = collect(1:n)
    irreversible_reactions_id, reversible_reactions_id = reversibility(lb, Reaction_Ids, printLevel)

    ## Create a new instance of the input model with homogenous bounds

    ModelObject_CC = Model_CC(S, Metabolites, Reactions, Genes, m, n, lb, ub)
    blocked_index, ν  = swiftCC(ModelObject_CC, SolverName, false, Tolerance, printLevel)
    blocked_index_rev = blocked_index ∩ reversible_reactions_id

    # Convert to Vector{Int64}:
    blocked_index_rev = convert(Vector{Int64}, blocked_index_rev)

    ## Correct Reversibility

    # Create an empty dictionary to store reaction bounds:
    Dict_bounds = Dict()

    # Iterate through each reaction to populate the dictionary with lower and upper bounds:
    for i = 1:n
        Dict_bounds[Reactions[i]] = lb[i], ub[i]
    end

    # Create a new Model_Correction object with the current data:
    ModelObject_Crrection = Model_Correction(S, Metabolites, Reactions, Genes, m, n, lb, ub, irreversible_reactions_id, reversible_reactions_id)

    # Apply distributedReversibility_Correction() to the model and update Reversibility, S and bounds:
    S, lb, ub, irreversible_reactions_id, reversible_reactions_id = distributedReversibility_Correction(ModelObject_Crrection, blocked_index_rev, SolverName, false)

    # Get the dimensions of the updated stoichiometric matrix:
    row, col = size(S)

    ## Count the number of reactions in each set

    n_irr = length(irreversible_reactions_id)
    n_rev = length(reversible_reactions_id)

    # Update the dictionary with the new lower and upper bounds after correction:
    for i = 1:n
        Dict_bounds[Reactions[i]] = lb[i], ub[i]
    end

    # Reconstruct the corrected model with updated parameters:
    model_Correction_Constructor(ModelObject_Crrection , S, Metabolites, Reactions, Genes, row, col, lb, ub, irreversible_reactions_id, reversible_reactions_id)

    ## Obtain blocked_index, fctable, Fc_Coefficients, and Dc_Coefficients

    ModelObject_QFCA = Model_QFCA(S, Metabolites, Reactions, Genes, row, col, lb, ub, irreversible_reactions_id, reversible_reactions_id)

    # Convert to Vector{Int64}:
    blocked_index = convert(Vector{Int64}, blocked_index)
    fctable, Fc_Coefficients, Dc_Coefficients = distributedQFCA(ModelObject_QFCA, blocked_index, SolverName, false)

    # Get the dimensions of fctable:
    row, col = size(fctable)

    # Remove blocked reactions from Reaction_Ids:
    Reaction_Ids_noBlocked = setdiff(Reaction_Ids, blocked_index)

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

    ## Iterate over sorted keys of FC_Final dictionary

    for key in sort(collect(keys(FC_Final)))
        # Append cluster members to FC_cluster_members array:
        append!(FC_cluster_members, FC_Final[key][2])

        # Append representative reaction to FC_representatives array:
        append!(FC_representatives, FC_Final[key][1])

        # Check if index_c is among the cluster members:
        if index_c in FC_Final[key][2]
            # If true, update index_c to the representative reaction:
            index_c = FC_Final[key][1]

            # Update Biomass variable with the reaction corresponding to index_c:
            Biomass = Reactions[index_c]
        end
    end

    ## Create a new dictionary FC_Clusters

    FC_Clusters = Dict()

    # Initialize counter variable:
    c = 1

    # Iterate over sorted keys of FC_Final dictionary:
    for key in sort(collect(keys(FC_Final)))
        # Assign each cluster from FC_Final to FC_Clusters with a numeric key:
        FC_Clusters[c] = FC_Final[key]

        # Increment counter for each cluster:
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

    # Sort the 'Reaction_Ids_noBlocked', 'irreversible_reactions_id', and 'reversible_reactions_id' arrays:
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

    # Check if we're using octuple precision (very high precision floating-point numbers):
    if OctuplePrecision
        # Define a model_irr using GenericModel from Clarabel.jl:
        model_local = GenericModel{BigFloat}(Clarabel.Optimizer{BigFloat})

        # Set verbose attribute to false (disable verbose output):
        set_attribute(model_local, "verbose", false)

        # Set absolute tolerance for gap convergence to 1e-32:
        set_attribute(model_local, "tol_gap_abs", 1e-32)

        # Set relative tolerance for gap convergence to 1e-32:
        set_attribute(model_local, "tol_gap_rel", 1e-32)
    else
        # If not using octuple precision, change the solver based on the solvername:
        model_local, solver = changeSparseQFCASolver(SolverName)
    end

    # Define the decision variables λ (for reactions), ν (for metabolites), and t (a scalar variable):
    @variable(model_local, λ[1:n])
    @variable(model_local, ν[1:m])
    @variable(model_local, t)

    # Set the objective function to minimize t (scalar variable):
    @objective(model_local, Min, t)

    # Constraint 1: λ and t must satisfy the NormOneCone (a type of norm constraint):
    con1 = @constraint(model_local, [t; λ] in MOI.NormOneCone(1 + length(λ)))

    # Constraint 2: The dual variable λ must be equal to the transposed stoichiometric matrix (S') times ν:
    con2 = @constraint(model_local, λ == S' * ν)

    # Constraint 3: Ensure that λ is 0 for reversible reactions (the reaction flux is zero for reversible reactions):
    con3 = @constraint(model_local, [j in reversible_reactions_id], λ[j] == 0.0)

    ## Perform distributed optimization for each reaction in the remove_list_DC list

    @sync @distributed for i in remove_list_DC

        # Save i as the row number to process the current reaction:
        row = i

        # Compute the linear combination for DCE by excluding the current reaction (i):
        DCE_LinearCombination = setdiff(1:n, i)

        ## Define additional constraints for the optimization problem

        # Constraint 4: Set λ to zero for specific reactions in the Eliminations ∩ DCE_LinearCombination set:
        con4 = @constraint(model_local, [j in Eliminations ∩ DCE_LinearCombination], λ[j] == 0.0)

        # Constraint 5: Set λ to be non-negative for the remaining reactions in DCE_LinearCombination:
        con5 = @constraint(model_local, [j in DCE_LinearCombination], λ[j] >= 0.0)

        # Constraint 6: Force λ for the current reaction (i) to be -1:
        con6 = @constraint(model_local, λ[i] == -1.0)

        # Solve the optimization problem for the current setup:
        optimize!(model_local)

        # Increment the counter to keep track of the number of optimizations performed:
        counter += 1

        # Get the values of the λ variables after solving the optimization problem:
        λ_vec = value.(λ)

        # Create an empty array to store the indices of the non-zero λ values:
        cols = Array{Int64}([])

        # Loop through the λ vector and check for values above the defined tolerance:
        for i = 1:length(λ_vec)
            if λ_vec[i] > Tolerance
                # Find the column index of the corresponding reaction in A_cols_reduced:
                index = findfirst(x -> x == i, A_cols_reduced)

                # Update the corresponding value in matrix A for the current row:
                A[row, index] = λ_vec[i]
            end
        end

        ## Condition 4: Remove all constraints in con4 from the model once used

        constraint_refs_con4 = [con4[i] for i in eachindex(con4)]
        for i in constraint_refs_con4
            delete(model_local, i)  # Remove the constraint from the model
            unregister(model_local, :i)  # Unregister the constraint for clean-up
        end

        ## Condition 5: Remove all constraints in con5 from the model once used

        constraint_refs_con5 = [con5[i] for i in eachindex(con5)]
        for i in constraint_refs_con5
            delete(model_local, i)  # Remove the constraint from the model
            unregister(model_local, :i)  # Unregister the constraint for clean-up
        end

        ## Condition 6: Remove the current λ constraint (con6) from the model

        delete(model_local, con6)  # Remove the specific λ[i] == -1 constraint
        unregister(model_local, :con6)  # Unregister this constraint for clean-up

    end

    # Convert matrix A to a matrix of Float64 data type for numerical stability:
    A = convert(Matrix{Float64}, A)

    # Get the number of rows and columns in matrix A:
    row_A, col_A = size(A)

    # Get the number of rows and columns in the stoichiometric matrix S:
    row_S, col_S = size(S)

    # Compute the modified stoichiometric matrix S̃ by multiplying S with A:
    S̃ = S * A

    # Remove any zero rows from the matrix S̃ and get the reduced metabolite list:
    S̃, Metabolites_reduced, Metabolites_elimination = remove_zeroRows(S̃, Metabolites)

    # Get the dimensions of the reduced stoichiometric matrix S̃:
    S̃_row, S̃_col = size(S̃)

    # Identify the reactions that have been eliminated (not part of A_cols_reduced):
    Reactions_elimination = Reactions[setdiff(range(1, n), A_cols_reduced)]

    # Get the reduced list of reactions corresponding to the columns in A_cols_reduced:
    R̃ = Reactions[A_cols_reduced]

    # Make a copy of the reduced metabolites to M̃:
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

    ## Genes

    # Iterate through each reaction in the model:
    for i in model.reactions
        # Iterate through the sorted reduction_map (which maps reactions to a list of associated reactions):
        for (key, value) in sort(reduction_map)
            # Check if the current reaction matches the first element in the reduction_map's value list:
            if i.first == Reactions[value[1]]
                # Loop through the second element of the value list, which contains a set of reactions:
                for j in value[2]
                    # Check if the gene association for the current reaction and the related reaction is not empty:
                    if !isnothing(i.second.gene_association_dnf) && !isnothing(model.reactions[Reactions[j]].gene_association_dnf)
                        # Concatenate the gene associations of the current reaction with those of the related reaction:
                        i.second.gene_association_dnf = vcat(i.second.gene_association_dnf, model.reactions[Reactions[j]].gene_association_dnf)
                    end
                end
            end
        end
    end

    # Iterate through each reaction in the model:
    for i in model.reactions
        # Check if the gene association exists and is not nothing:
        if !isnothing(i.second.gene_association_dnf)
            # Flatten the non-empty elements of the gene_association_dnf list into a single iterable:
            flattened = filter(!isempty, i.second.gene_association_dnf) |> Iterators.flatten
            # Create a set to keep track of unique genes:
            seen = Set{String}()
            # Create an empty vector to store unique gene associations:
            result = Vector{Vector{String}}()
            # Iterate through the flattened list of genes:
            for gene in flattened
                # Check if the gene has not been seen before:
                if !(gene in seen)
                    # Add the gene to the set of seen genes:
                    push!(seen, gene)
                    # Add the gene as a new vector to the result list:
                    push!(result, [gene])
                end
            end
            # Update the gene_association_dnf with the unique gene associations:
            i.second.gene_association_dnf = result
        end
    end

    ## Update lb & Up

    # Iterate through each reaction in the model to update the lower and upper bounds:
    for i in model.reactions
        # Set the lower bound for the reaction using the pre-calculated Dict_bounds:
        i.second.lower_bound = Dict_bounds[i.first][1]
        # Set the upper bound for the reaction using the pre-calculated Dict_bounds:
        i.second.upper_bound = Dict_bounds[i.first][2]
    end

    # Filter out reactions that are in the list of eliminated reactions:
    filter!(pair -> !(pair.first in Reactions_elimination), model.reactions)

    # Filter out metabolites that are in the list of eliminated metabolites:
    filter!(pair -> !(pair.first in Metabolites_elimination), model.metabolites)

    ## Update Stoichiometry Matrix

    # Iterate through each reaction in the model to clear the stoichiometry:
    for i in model.reactions
        # Collect the keys from the stoichiometry map of the current reaction:
        for key in collect(keys(i.second.stoichiometry))
            # Remove the stoichiometry entry for the current key:
            delete!(i.second.stoichiometry, key)
        end
    end

    # Iterate through each reaction to update the stoichiometry based on the modified matrix:
    for i in model.reactions
        # Find the index of the current reaction in the reduced reaction list R̃:
        index_col = findfirst(x -> x == i.first, R̃)

        # Get the corresponding column in the modified stoichiometric matrix S̃:
        stoichiometry_vector = S̃[:,index_col]

        # Initialize a counter for metabolites:
        met = 1

        # Loop through the stoichiometry vector to update the stoichiometry for each metabolite:
        for c = 1:length(stoichiometry_vector)
            # Push the metabolite and its corresponding stoichiometry value to the reaction:
            push!(i.second.stoichiometry, "$(M̃[met])" => stoichiometry_vector[c])
            met += 1
        end
    end

    ## Update Biomass

    # Iterate through the reduction map to check and update the Biomass reaction:
    for (key, value) in sort(reduction_map)
        # If the biomass reaction (index_c) is in the reduction map, update its index:
        if index_c in reduction_map[key][2]
            index_c = reduction_map[key][1]
            Biomass = Reactions[index_c]
        end
    end

    ## Update Objective coefficient

    # Iterate through the reactions to set the objective coefficient for the biomass reaction:
    for i in model.reactions
        # If the current reaction matches the biomass reaction, set its objective coefficient to 1.0:
        if i.first == Biomass
            i.second.objective_coefficient = 1.0
        end
    end

    ## Update Genes

    # Initialize an empty list to store final gene associations:
    Genes_final = []

    # Iterate through each reaction to collect gene associations:
    for i in model.reactions
        # If the reaction has a non-empty gene association, append it to Genes_final:
        if !isnothing(i.second.gene_association_dnf)
            append!(Genes_final, i.second.gene_association_dnf)
        end
    end

    # Remove duplicate genes from Genes_final:
    Genes_final = unique(Genes_final)

    # Flatten the nested gene lists in Genes_final into a single list:
    Genes_final = collect(Iterators.flatten(Genes_final))

    # Determine the genes that need to be removed by finding the difference:
    Genes_removal = setdiff(Genes, Genes_final)

    # Filter out genes in the model that are present in Genes_removal:
    filter!(pair -> !(pair.first in Genes_removal), model.genes)

    ## Print out results if requested

    if printLevel > 0
        printstyled("Metabolic Network Reductions:\n"; color=:cyan)
        if OctuplePrecision
            printstyled("Solver = Clarabel \n"; color=:green)
        else
            printstyled("Solver = $SolverName\n"; color=:green)
        end
        printstyled("Tolerance = $Tolerance\n"; color=:magenta)
        println("Original Network:")
        println("S           : $(row_S) x $(col_S)")
        println("Genes       : $(length(Genes))")
        println("Metabolites : $(m)")
        println("Reactions   : $(n)")
        println("Reduced Network:")
        println("S           : $(S̃_row) x $(S̃_col)")
        println("Genes       : $(length(Genes_final))")
        println("Metabolites : $(length(Metabolites_reduced))")
        println("Reactions   : $(length(R̃))")
        println("A matrix    : $(row_A) x $(col_A)")
    end
    return model
end

end

#-----------------------------------------------------------------------------------------------------------------------------------------------------
