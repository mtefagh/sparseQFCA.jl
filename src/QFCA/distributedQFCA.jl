#-----------------------------------------------------------------------------------------------------------------------------------------------------
#=
    Purpose:    Parallel computation of quantitative flux coupling using swiftCC algorithm and Gaussian Elimination
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology
    Date:       October 2022
=#
#-----------------------------------------------------------------------------------------------------------------------------------------------------

module DistributedQFCA

export Model_QFCA, model_QFCA_Constructor, distributedQFCA

using JuMP, COBREXA, LinearAlgebra, SparseArrays, Distributed, SharedArrays

include("../Pre_Processing/Pre_processing.jl")

using .Pre_processing

include("../ConsistencyChecking/SwiftCC.jl")

using .SwiftCC

"""
    Model_QFCA(S, Metabolites, Reactions, Genes, m, n, lb, ub, irreversible_reactions_id, reversible_reactions_id)

A general type for storing a CanonicalModel which contains the following fields to run distributedQFCA():

-`ModelObject_QFCA`:                A newly object of Model_Correction.
- `S`:                              LHS matrix (m x n)
- `Metabolites`:                    List of metabolic network metabolites.
- `Reactions`:                      List of metabolic network reactions.
- `Genes`:                          List of metabolic network reactions.
- `m`:                              Number of rows of stoichiometric matrix.
- `n`:                              Number of columns of stoichiometric matrix.
- `lb`:                             Lower bound vector (n x 1)
- `ub`:                             Upper bound vector (n x 1)
- `irreversible_reactions_id`:      IDs of irreversible reactions.
- `reversible_reactions_id`  :      IDs of reversible reactions.

"""

#-------------------------------------------------------------------------------------------

mutable struct Model_QFCA
    S                               ::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}
    Metabolites                     ::Array{String,1}
    Reactions                       ::Array{String,1}
    Genes                           ::Array{String,1}
    m                               ::Int
    n                               ::Int
    lb                              ::Array{Float64,1}
    ub                              ::Array{Float64,1}
    irreversible_reactions_id       ::Array{Int64}
    reversible_reactions_id         ::Array{Int64}
end

#-------------------------------------------------------------------------------------------

"""
    model_QFCA_Constructor(ModelObject, S, Metabolites, Reactions, Genes, m, n, lb, ub, irreversible_reactions_id, reversible_reactions_id)

The function takes in several arguments, including a ModelObject of type Model_QFCA,
and assigns values to its fields based on the other arguments passed in.

# INPUTS

-`ModelObject_QFCA`:                A newly object of Model_Correction.
- `S`:                              LHS matrix (m x n)
- `Metabolites`:                    List of metabolic network metabolites.
- `Reactions`:                      List of metabolic network reactions.
- `Genes`:                          List of metabolic network reactions.
- `m`:                              Number of rows of stoichiometric matrix.
- `n`:                              Number of columns of stoichiometric matrix.
- `lb`:                             Lower bound vector (n x 1)
- `ub`:                             Upper bound vector (n x 1)
- `irreversible_reactions_id`:      IDs of irreversible reactions.
- `reversible_reactions_id`  :      IDs of reversible reactions.

"""

function model_QFCA_Constructor(ModelObject::Model_QFCA, S::Union{SparseMatrixCSC{Float64,Int64}, AbstractMatrix}, Metabolites::Array{String,1}, Reactions::Array{String,1},
                                Genes::Array{String,1},m::Int, n::Int, lb::Array{Float64,1}, ub::Array{Float64,1}, irreversible_reactions_id::Array{Int64}, reversible_reactions_id::Array{Int64})
     ModelObject_QFCA.S = S
     ModelObject_QFCA.Metabolites = Metabolites
     ModelObject_QFCA.Reactions = Reactions
     ModelObject_QFCA.Genes = Genes
     ModelObject_QFCA.m = m
     ModelObject_QFCA.n = n
     ModelObject_QFCA.lb = lb
     ModelObject_QFCA.ub = ub
     ModelObject_QFCA.irreversible_reactions_id = irreversible_reactions_id
     ModelObject_QFCA.reversible_reactions_id = reversible_reactions_id
end


"""
    distributedQFCA(myModel)

The function first exports data from ModelObject_QFCA and calculates directional couplings and coefficients for each pair of reactions in the model.
It also has the option to remove reactions from the model and recalculate the couplings. The function uses parallel processing to distribute
the calculations across multiple processors. The output is a matrix that shows the coupling relations between the reactions in the model.

# INPUTS

-`ModelObject_QFCA`:           A newly object of Model_QFCA.
- `blocked_index`:             IDs of of blocked reactions.

# OPTIONAL INPUTS

- `SolverName`:                Name of the solver(default: HiGHS).
- `OctuplePrecision`:          A flag(default: false) indicating whether octuple precision should be used when solving linear programs.
- `removing`:                  A boolean variable that indicates whether reactions should be removed from the network in the stages of determining coupling or not.
- `Tolerance`:                 A small number that represents the level of error tolerance.
- `printLevel`:                Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

- `fctable`:                   The resulting flux coupling matrix.
                               The meaning of the entry (i, j) is:
                                    * 0 - uncoupled reactions
                                    * 1 - fully coupled reactions
                                    * 2 - partially coupled reactions
                                    * 3 - reaction i is directionally coupled to reaction j
                                    * 4 - reaction j is directionally coupled to reaction i
- `Fc_Coefficients`:           A list of fully-coupling coefficients for each reaction in the model.
- `Dc_Coefficients`:           A list of DCE (directional coupling equation) coefficients for each reaction in the model.

# EXAMPLES

- Full input/output example
```julia
julia> fctable, Fc_Coefficients, Dc_Coefficients = distributedQFCA(ModelObject_QFCA, blocked_index)
```

See also: `dataOfModel()`, `reversibility()`, `homogenization()`, `Model_QFCA`, `model_QFCA_Constructor()`, `distributedReversibility_Correction()`

"""

function distributedQFCA(ModelObject_QFCA::Model_QFCA, blocked_index::Vector{Int64}, SolverName::String="HiGHS", OctuplePrecision::Bool=false, removing::Bool=false, Tolerance::Float64=1e-6, printLevel::Int=1)

    ## Extract relevant information from the ModelObject_QFCA

    # Extracting the stoichiometric matrix (S) from the ModelObject_QFCA:
    S = ModelObject_QFCA.S

    # Extracting the array of metabolite IDs from the ModelObject_QFCA:
    Metabolites = ModelObject_QFCA.Metabolites

    # Extracting the array of reaction IDs from the ModelObject_QFCA:
    Reactions = ModelObject_QFCA.Reactions

    # Extracting the array of gene IDs from the ModelObject_QFCA:
    Genes = ModelObject_QFCA.Genes

    # Extracting the number of metabolites (m) from the ModelObject_QFCA:
    m = ModelObject_QFCA.m

    # Extracting the number of reactions (n) from the ModelObject_QFCA:
    n = ModelObject_QFCA.n

    # Extracting the lower bounds for reactions from the ModelObject_QFCA:
    lb = ModelObject_QFCA.lb

    # Extracting the upper bounds for reactions from the ModelObject_QFCA:
    ub = ModelObject_QFCA.ub

    # Extracting the IDs of irreversible reactions from the ModelObject_QFCA:
    irreversible_reactions_id = ModelObject_QFCA.irreversible_reactions_id

    # Extracting the IDs of reversible reactions from the ModelObject_QFCA:
    reversible_reactions_id = ModelObject_QFCA.reversible_reactions_id

    ## Count the number of reactions in each set

    n_irr = length(irreversible_reactions_id)
    n_rev = length(reversible_reactions_id)

    # Remove blocked reactions from Reaction_Ids:
    Reaction_Ids = collect(1:n)
    Reaction_Ids_noBlocked = setdiff(Reaction_Ids, blocked_index)

    # Calculate the number of blocked reactions:
    n_blocked = length(blocked_index)

    # Create a new array of reaction IDs that excludes the blocked reactions:
    Reactions_noBlocked = Reactions[setdiff(range(1, n), blocked_index)]

    # Create a new array of lower bounds for the non-blocked reactions:
    lb_noBlocked = lb[setdiff(1:end, blocked_index)]

    # Create a new array of upper bounds for the non-blocked reactions:
    ub_noBlocked = ub[setdiff(1:end, blocked_index)]

    # Create a submatrix of the stoichiometric matrix 'S' with only non-blocked reactions:
    S_noBlocked = S[:, setdiff(1:end, blocked_index)]

    # Update the array of irreversible reaction IDs by removing the blocked reactions:
    irreversible_reactions_id = setdiff(irreversible_reactions_id, blocked_index)

    # Update the array of reversible reaction IDs by removing the blocked reactions:
    reversible_reactions_id = setdiff(reversible_reactions_id, blocked_index)

    # Calculate the new number of irreversible reactions after excluding blocked ones:
    n_irr = length(irreversible_reactions_id)

    # Calculate the new number of reversible reactions after excluding blocked ones:
    n_rev = length(reversible_reactions_id)

    # Convert to Vector{Int64}:
    irreversible_reactions_id = convert(Vector{Int64}, irreversible_reactions_id)
    # Convert to Vector{Int64}:
    reversible_reactions_id = convert(Vector{Int64}, reversible_reactions_id)

    # Get the number of rows(metabolites) and columns(reactions) in the S_noBlocked matrix:
    row_noBlocked, col_noBlocked = size(S_noBlocked)

    # Create a new array of lower bounds for the non-blocked reactions:
    lb_noBlocked = lb[setdiff(1:end, blocked_index)]

    # Create a new array of upper bounds for the non-blocked reactions:
    ub_noBlocked = ub[setdiff(1:end, blocked_index)]

    ## Remove all rows from a given sparse matrix S that contain only zeros and the corresponding metabolites from the Metabolites array

    S_noBlocked, Metabolites = remove_zeroRows(S_noBlocked,Metabolites)

    ## Create a new instance of the input model with homogenous bounds

    ModelObject_CC = Model_CC(S_noBlocked, Metabolites, Reactions_noBlocked, Genes, row_noBlocked, col_noBlocked, lb_noBlocked, ub_noBlocked)

    ## Create an empty matrix to store directional couplings

    DC_Matrix = SharedArray{Int,2}((col_noBlocked, col_noBlocked), init = false)

    ## Create an empty matrix to store directional coupling coefficients

    Dc_Coefficients = SharedArray{Float64,2}((col_noBlocked, col_noBlocked), init = false)

    ## Calculate directional couplings and coefficients for each pair of reactions

    # This is a parallel loop that distributes iterations across multiple processors:
    @sync @distributed for i = 1:col_noBlocked
        if(removing)
            # Store current values of bounds and stoichiometric matrix before removing the ith reaction from the metabolic network:
            lb_noBlocked_temp = copy(lb_noBlocked)
            ub_noBlocked_temp = copy(ub_noBlocked)
            S_noBlocked_temp  = copy(S_noBlocked)

            # Remove the ith reaction from the metabolic network:
            lb_noBlocked = lb_noBlocked[setdiff(1:end, i)]
            ub_noBlocked = ub_noBlocked[setdiff(1:end, i)]
            S_noBlocked  = S_noBlocked[:, setdiff(1:end, i)]
            row_noBlocked, col_noBlocked = size(S_noBlocked)

            # Calculate the set of blocked reactions and dual variables for the modified network:
            model_CC_Constructor(ModelObject_CC ,S_noBlocked, Metabolites, Reactions_noBlocked, Genes, row_noBlocked, col_noBlocked, lb_noBlocked, ub_noBlocked)
            blocked, ν, status = swiftCC(ModelObject_CC, SolverName, false, Tolerance, 0)

            # Update indices of blocked reactions after removing ith reaction:
            for j = 1:length(blocked)
                if blocked[j] >= i
                    blocked[j] = blocked[j] + 1
                end
            end

            # Update DC_Matrix based on the blocked reactions:
            DC_Matrix[i,blocked] .= 1.0

            # Calculate the coefficients for DC(i, j) for all blocked reactions j:
            λ = S_noBlocked_temp' * ν
            Dc_Coefficients[i, :] = λ

            # Retrieve the original values of the bounds and stoichiometric matrix:
            lb_noBlocked = copy(lb_noBlocked_temp)
            ub_noBlocked = copy(ub_noBlocked_temp)
            S_noBlocked  = copy(S_noBlocked_temp)
        else
            # Store current values of bounds before resetting the bounds of the ith reaction to zero:
            lb_temp = lb_noBlocked[i]
            ub_temp = ub_noBlocked[i]

            # Reset the bounds of the ith reaction to zero:
            lb_noBlocked[i] = 0.0
            ub_noBlocked[i] = 0.0

            # Find the set of blocked reactions and dual variables for the modified network:
            model_CC_Constructor(ModelObject_CC ,S_noBlocked, Metabolites, Reactions_noBlocked, Genes, row_noBlocked, col_noBlocked, lb_noBlocked, ub_noBlocked)
            blocked, ν = swiftCC(ModelObject_CC, SolverName, false, Tolerance, 0)

            # Update DC_Matrix based on the blocked reactions:
            DC_Matrix[i,blocked] .= 1.0

            # Calculate the coefficients for DC(i, j) for all blocked reactions j:
            λ = S_noBlocked' * ν
            Dc_Coefficients[i, :] = λ

            # Retrieve the original values of the bounds:
            lb_noBlocked[i] = lb_temp
            ub_noBlocked[i] = ub_temp
        end
    end

    ## Define a matrix to save all coupling relations

    fctable = SharedArray{Int,2}((col_noBlocked, col_noBlocked), init = false)

    ## Directional Coupling

    ## Update the functional coupling table based on the directional coupling matrix

    # For each reaction i and reaction j in the range of blocked reactions:
    for i = 1:col_noBlocked
        for j = 1:col_noBlocked
            # If there is directional coupling from reaction j to reaction i:
            if DC_Matrix[i, j] == 1.0
                fctable[j,i] = 3.0
            end
        end
    end

    ## Self-Coupling

    fctable[diagind(fctable)] .= 1.0

    ## Distinguish between (DC, FC or PC) kinds of coupling

    for i = 1:col_noBlocked
        for j = i+1:col_noBlocked
            # If i has an arrow pointing to j, but j doesn't have an arrow pointing to i, add an arrow from j to i:
            if fctable[i,j] == 3.0 && fctable[j,i] == 0.0
                fctable[j,i] = 4.0

            # If j has an arrow pointing to i, but i doesn't have an arrow pointing to j, add an arrow from i to j:
            elseif fctable[i,j] == 0.0 && fctable[j,i] == 3.0
                fctable[i,j] = 4.0

            # If i and j have arrows pointing in both directions, set the arrows to be bidirectional:
            elseif fctable[i,j] == 3.0 && fctable[j,i] == 3.0
                fctable[i,j] = fctable[j,i] = 2.0

            # If there are no arrows between i and j, or both arrows are already bidirectional, do nothing and continue to the next pair of reactions:
            else
                continue
            end
        end
    end

    # Create a dictionary to store partially coupled pairs:
    PC = Dict()

    # Initialize a variable to keep track of the partial couples:
    s = 1

    # Initialize an empty array to store partially coupled pairs:
    Partially_Couples = []

    ## Iterate over all pairs of nodes in the coupling table

    for i = 1:col_noBlocked
        for j = 1:col_noBlocked
            # If the nodes are partially coupled in both directions:
            if fctable[i,j] == fctable[j,i] == 2.0

                # Save the partially coupled pair:
                Partially_Couples = i,j

                # Add the pair to the dictionary of partially coupled pairs:
                PC[s] = Partially_Couples

                # Increment the counter for partially coupled pairs:
                s = s + 1
            end
        end
    end

    # Create a matrix to store the fully coupled coefficients:
    Fc_Coefficients = SharedArray{Float64,2}((col_noBlocked, col_noBlocked), init = false)

    # Set the self-coupling coefficients to 1:
    Fc_Coefficients[diagind(Fc_Coefficients)] .= 1.0

    ## Create a dictionary to store fully coupled by one Metabolite concept

    FC_OneMet = Dict()

    # Initialize a variable to keep track of fully coupled by one Metabolite:
    met = 1

    ## FC_OneMet    Key : #Met  ---> Values : Indices of non-zero values

    ## Loop through each row in S matrix

    for row ∈ eachrow(S_noBlocked)
        non_zero_indices = []
        row = sparsevec(row)

        # Check if row has exactly 2 non-zero values:
        if nnz(row) == 2
            # Loop through row and get the indices of the non-zero values:
            for i = 1:length(row)
                if row[i] != 0
                    append!(non_zero_indices, i)
                end
            end

            # Add the non-zero indices as a tuple to FC_OneMet dictionary with the key met:
            FC_OneMet[met] = (non_zero_indices[1], non_zero_indices[2])
        end

        # Increment met counter:
        met += 1
    end

    ## Convert FC_OneMet : Vector{any} ---> Tuple{Int64, Int64}

    # Loop through each key in FC_OneMet dictionary and convert the values from vector to tuple:
    for key ∈ sort(collect(keys(FC_OneMet)))
        FC_OneMet[key] = (FC_OneMet[key][1], FC_OneMet[key][2])
        FC_OneMet[key] = convert(Tuple{Int64, Int64}, FC_OneMet[key])
    end

    ## Loop through each key in FC_OneMet dictionary

    for key ∈ sort(collect(keys(FC_OneMet)))
        # Check if the row index and column index of the key in PC matches the tuple value in FC_OneMet:
        if (PC[key][1] == FC_OneMet[key][1]) && (PC[key][2] == FC_OneMet[key][2])

        # Edit fctable and changing partially Coupling to Fully Coupling
        # Set the corresponding entries in the coupling table to 1:
        fctable[PC[key][1],PC[key][2]] = 1.0
        fctable[PC[key][2],PC[key][1]] = 1.0

        # Calculate the fully coupled coefficients and set the values in Fc_Coefficients matrix:
        Fc_Coefficients[PC[key][1],PC[key][2]] = abs(S_noBlocked[key ,FC_OneMet[key][2]]) / abs(S_noBlocked[key ,FC_OneMet[key][1]])
        Fc_Coefficients[PC[key][2],PC[key][1]] = abs(S_noBlocked[key ,FC_OneMet[key][1]]) / abs(S_noBlocked[key ,FC_OneMet[key][2]])
        end
    end

    ## Solve an LU decomposition for each pair of partially coupled nodes to determine fully coupling

    # Start a distributed loop over keys in PC, sorted in ascending order:
    @sync @distributed for key ∈ sort(collect(keys(PC)))

        ## Check if the value of fctable at the location given by the key in PC is not equal to 1.0

        if fctable[PC[key][1],PC[key][2]] != 1.0

            # Transpose S_noBlocked:
            S_noBlocked_transpose = S_noBlocked'

            # Initialize CO as a zeros vector with length equal to col_noBlocked:
            CO = zeros(col_noBlocked)

            # Set the value of CO at the index given by the first element of the value of key in PC to -1:
            CO[PC[key][1]] = -1

            # Remove the row given by the second element of the value of key in PC from S_noBlocked_transpose:
            S_noBlocked_transpose_removedC = S_noBlocked_transpose[1:end .!= PC[key][2], :]

            # Remove the element at the index given by the second element of the value of key in PC from CO:
            CO = CO[setdiff(1:end, PC[key][2])]

            ## Use Gaussian elimination to solve for X in the equation S_noBlocked_transpose_removedC * X = CO

            X = S_noBlocked_transpose_removedC \ CO

            ## Calculate the solution Sol by subtracting CO from S_noBlocked_transpose_removedC * X

            Sol = (S_noBlocked_transpose_removedC*X) - (CO)

            ## Check if the solution vector "Sol" is close to the zero vector with a tolerance of "Tolerance"

            if isapprox(norm(Sol), 0.0, atol = Tolerance)

                # Edit fctable and change partially Coupling to Fully Coupling:
                fctable[PC[key][1],PC[key][2]] = 1.0
                fctable[PC[key][2],PC[key][1]] = 1.0

                # Cast SparseVector{Float64, Int64} to SparseMatrixCSC{Float64, Int64}:
                rowC_S_noBlocked_transpose = sparsevec(S_noBlocked_transpose[PC[key][2], :])
                matrix_rowC_S_noBlocked_transpose = sparse(rowC_S_noBlocked_transpose)
                sparseMatrix_rowC = SparseMatrixCSC(matrix_rowC_S_noBlocked_transpose)

                # Calculate coefficient:
                C = sparseMatrix_rowC' * X

                ## Set coefficient for i,j

                # V_i = C * V_j:
                Fc_Coefficients[PC[key][1],PC[key][2]] = C[1]

                # V_j = 1/C * V_i:
                Fc_Coefficients[PC[key][2],PC[key][1]] = 1 / C[1]
            end
        end
    end

    ## Print out results if requested

    if printLevel > 0
        d_0 = 0
        d_1 = 0
        d_2 = 0
        d_3 = 0
        d_4 = 0
        d_0 = sum(fctable .== 0.0)
        d_1 = sum(fctable .== 1.0)
        d_2 = sum(fctable .== 2.0)
        d_3 = sum(fctable .== 3.0)
        d_4 = sum(fctable .== 4.0)
        printstyled("Distributed Quantitative Flux Coupling Analysis(distributedQFCA):\n"; color=:cyan)
        println("Number of Proccess : $(nprocs())")
        println("Number of Workers  : $(nworkers())")
        if OctuplePrecision
            printstyled("The name of the solver = Clarabel \n"; color=:green)
        else
            printstyled("The name of the solver = $SolverName\n"; color=:green)
        end
        printstyled("Tolerance = $Tolerance\n"; color=:magenta)
        println("Final fctable : ")
        println("Number of 0's (unCoupled) : $d_0")
        println("Number of 1's (Fully)     : $d_1")
        println("Number of 2's (Partialy)  : $d_2")
        println("Number of 3's (DC i-->j)  : $d_3")
        println("Number of 4's (DC j-->i)  : $d_4")
    end

    # Convert fctable to a matrix of integers:
    fctable = convert(Matrix{Int}, fctable)

    # Convert Fc_Coefficients to a matrix of Float64:
    Fc_Coefficients = convert(Matrix{Float64}, Fc_Coefficients)

    # Convert Dc_Coefficients to a matrix of Float64:
    Dc_Coefficients = convert(Matrix{Float64}, Dc_Coefficients)

    return fctable, Fc_Coefficients, Dc_Coefficients
    end

end

#-----------------------------------------------------------------------------------------------------------------------------------------------------
