#-------------------------------------------------------------------------------------------
#=
    Purpose:    QFCA computes the table of flux coupling relations for a metabolic network
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology - Iran
    Date:       October 2022
=#
#-------------------------------------------------------------------------------------------

module DistributedQFCA
export addQFCAProcs, removeQFCAProcs, distributedQFCA

using Distributed

@everywhere using GLPK, JuMP, COBREXA, LinearAlgebra, SparseArrays, Distributed, SharedArrays

include("../Data Processing/pre_processing.jl")

using .pre_processing

include("../Consistency Checking/SwiftCC.jl")

using .SwiftCC

"""
    addQFCAProcs(n)

A function that adds n-1 processes to System.

# INPUTS

- `n`:        Number of processes to add.

# OPTIONAL INPUTS

-

# OUTPUTS

-

# EXAMPLES

- Full input/output example
```julia
julia> addQFCAProcs(n)
```

See also: ``

"""

function addQFCAProcs(n::Int)
   addprocs(n-1)
   return
end

"""
    removeQFCAProcs()

A function that removes n-1 processes from Systems.

# INPUTS

-

# OPTIONAL INPUTS

-

# OUTPUTS

-

# EXAMPLES

- Full input/output example
```julia
julia> removeQFCAProcs()
```

See also: ``

"""

function removeQFCAProcs()
    procs_vector = procs()
    for i in procs_vector[2:end]
       rmprocs(i)
    end
    return
end

"""
    distributedQFCA(myModel)

A function that computes the table of flux coupling relations for a metabolic network.

# INPUTS

- `myModel`:        A model that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

-

# OUTPUTS

- `fctable`:                   The resulting flux coupling matrix.
                               The meaning of the entry (i, j) is:
                                    * 0 - uncoupled reactions
                                    * 1 - fully coupled reactions
                                    * 2 - partially coupled reactions
                                    * 3 - reaction i is directionally coupled to reaction j
                                    * 4 - reaction j is directionally coupled to reaction i
- `Fc_Coefficients`:           fully-coupling coefficients
- `Dc_Coefficients`:           DCE coefficients.
- `n_blocked`:                 Number of blocked reactions.

# EXAMPLES

- Full input/output example
```julia
julia> fctable = distributedQFCA(myModel)
```

See also: `dataOfModel()`, `reversibility()`, `homogenization()`, `MyModel`, `myModel_Constructor()`

"""

function distributedQFCA(myModel::StandardModel, removing::Bool=false)

    # Exporting data from StandardModel:

    S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel)

    # Determining the reversibility of a reaction:

    irreversible_reactions_id, reversible_reactions_id = reversibility(lb)

    # Determining the number of irreversible and reversible reactions:

    n_irr = length(irreversible_reactions_id)
    n_rev = length(reversible_reactions_id)

    # Homogenizing the upper_bound and lower_bound of reactions:

    lb, ub = homogenization(lb, ub)

    # assigning a small value to Tolerance representing the level of error tolerance:

    Tolerance = getTolerance()

    # Creating a newly object of MyModel:

    ModelObject = MyModel(S, Metabolites, Reactions, Genes, m, n, lb, ub)

    # Removing blocked reactions from metabolic network:

    blocked_index, dualVar  = swiftCC(ModelObject)
    n_blocked = length(blocked_index)
    Reactions_noBlocked = Reactions[setdiff(range(1, n), blocked_index)]
    lb_noBlocked = lb[setdiff(1:end, blocked_index)]
    ub_noBlocked = ub[setdiff(1:end, blocked_index)]
    S_noBlocked  = S[:, setdiff(1:end, blocked_index)]
    row_noBlocked, col_noBlocked = size(S_noBlocked)

    # Determining Coupling by using swiftCC function:

    @everywhere D = Dict()

    # Defining a matrix to save DC-couplings:

    DC_Matrix = SharedArray{Int,2}((col_noBlocked, col_noBlocked), init = false)

    # Defining a matrix to save DCE coefficients:

    Dc_Coefficients = SharedArray{Int,2}((col_noBlocked, col_noBlocked), init = false)

    # Finding Coupling between reactions by using swiftCC functions in Parallel:

    @sync @distributed for i in range(1, col_noBlocked)

        if(removing)

            # Saving Values:

            lb_noBlocked_temp = copy(lb_noBlocked)
            ub_noBlocked_temp = copy(ub_noBlocked)
            S_noBlocked_temp = copy(S_noBlocked)

            # Removing:

            lb_noBlocked = lb_noBlocked[setdiff(1:end, i)]
            ub_noBlocked = ub_noBlocked[setdiff(1:end, i)]
            S_noBlocked  = S_noBlocked[:, setdiff(1:end, i)]
            row_noBlocked, col_noBlocked = size(S_noBlocked)

            # Finding Couples:

            myModel_Constructor(ModelObject ,S_noBlocked, Metabolites, Reactions_noBlocked, Genes, row_noBlocked, col_noBlocked, lb_noBlocked, ub_noBlocked)
            blocked, dualVar = swiftCC(ModelObject, epsilon)
            for j = 1:length(blocked)
                if blocked[j] >= i
                   blocked[j] = blocked[j] + 1
                end
            end

            D[i] = blocked
            D_values = collect(values(D[i]))
            DC_Matrix[i,D_values] .= 1.0

            lambda = S_noBlocked_temp' * dualVar
            Dc_Coefficients[i, :] = lambda

            # Retrieving Values

            lb_noBlocked = copy(lb_noBlocked_temp)
            ub_noBlocked = copy(ub_noBlocked_temp)
            S_noBlocked = copy(S_noBlocked_temp)
        else

            # Saving Values

            lb_temp = lb_noBlocked[i]
            ub_temp = ub_noBlocked[i]

            # Reset to zero

            lb_noBlocked[i] = 0.0
            ub_noBlocked[i] = 0.0

            # Finding Couples

            myModel_Constructor(ModelObject ,S_noBlocked, Metabolites, Reactions_noBlocked, Genes, row_noBlocked, col_noBlocked, lb_noBlocked, ub_noBlocked)
            blocked, dualVar = swiftCC(ModelObject, epsilon)

            D[i] = blocked
            D_values = collect(values(D[i]))
            DC_Matrix[i,D_values] .= 1.0

            lambda = S_noBlocked' * dualVar
            Dc_Coefficients[i, :] = lambda

            # Retrieving Values

            lb_noBlocked[i] = lb_temp
            ub_noBlocked[i] = ub_temp
        end
    end

    # Defining a matrix to save all coupling relations:

    fctable = SharedArray{Int,2}((col_noBlocked, col_noBlocked), init = false)

    # Directional Coupling:

    # i-->j

    for i in range(1, col_noBlocked)
        for j in range(1, col_noBlocked)
            if DC_Matrix[i, j] == 1
                fctable[j,i] = 3.0
            end
        end
    end

    # Self-Coupling:

    fctable[diagind(fctable)] .= 1.0

    # Distingushing between (DC, FC or PC) kinds of coupling:

     for i in range(1, col_noBlocked)
        for j in range(i+1, col_noBlocked)

            # i--->j

            if fctable[i,j] == 3.0 && fctable[j,i] == 0.0
                fctable[j,i] = 4.0

            # j--->i

           elseif fctable[i,j] == 0.0 && fctable[j,i] == 3.0
                fctable[i,j] = 4.0

            # i--->j && j--->i   ----> Partially Coupling

            elseif fctable[i,j] == 3.0 && fctable[j,i] == 3.0
                fctable[i,j] = fctable[j,i] = 2.0
            else
                continue
            end
        end
    end

    # Determining Partial Couples:

    @everywhere PC = Dict()

    s = 1
    Partially_Couples = []

    for i in range(1, col_noBlocked)
        for j in range(1, col_noBlocked)
            if fctable[i,j] == fctable[j,i] == 2.0
                Partially_Couples = i,j
                PC[s] = Partially_Couples
                s = s + 1
            end
        end
    end

    # Defining a matrix to save fully-coupling coefficients:

    Fc_Coefficients = SharedArray{Float64,2}((col_noBlocked, col_noBlocked), init = false)

    # Self-Coupling coefficients:
    #   V_i = 1 * V_i

    Fc_Coefficients[diagind(Fc_Coefficients)] .= 1.0

    # Defining a dictionary to save metabolites with only two reactions that are fully-coupled:
    @everywhere FC_OneMet = Dict()
    met = 1
    for row in eachrow(S)
        non_zero_indices = []
        row = sparsevec(row)
        if nnz(row) == 2
            for i = 1:length(row)
                if row[i] != 0
                append!(non_zero_indices, i)
                end
            end
            FC_OneMet[met] = non_zero_indices
        end
        met += 1
    end

    # Converting FC_OneMet : Vector{any} ---> Tuple{Int64, Int64}

    for key in sort(collect(keys(FC_OneMet)))
        FC_OneMet[key] = (FC_OneMet[key][1], FC_OneMet[key][2])
        FC_OneMet[key] = convert(Tuple{Int64, Int64}, FC_OneMet[key])
    end

    # Merging Keys of two Dictionary to use it as a Tuple:

    Keys = Tuple(zip(sort(collect(keys(PC))), sort(collect(keys(FC_OneMet)))))

    # Solving a LU for each pair to determine Fully Coupling:

    @sync @distributed for (key_PC, Key_FC_OneMet) in Keys

        # Check if the nodes corresponding to the current key_PC in the Partially-Coupling Dictionary (PC) and
        # the nodes corresponding to the current Key_FC_OneMet in the Fully Connected One-Metabolite Dictionary (FC_OneMet)
        # are the same.

        if (PC[key_PC][1] == FC_OneMet[Key_FC_OneMet][1]) && (PC[key_PC][2] == FC_OneMet[Key_FC_OneMet][2])

            # If the nodes are the same, update the FC table to indicate that the nodes are fully coupled.

            fctable[PC[key_PC][1],PC[key_PC][2]] = 1.0
            fctable[PC[key_PC][2],PC[key_PC][1]] = 1.0

            # Calculate the fully coupling coefficients using the ratio of the absolute values of the
            # corresponding elements in the No-Blocked Stoichiometric matrix (S_noBlocked).

            # V_i = |b/a| * V_j

            Fc_Coefficients[PC[key_PC][1],PC[key_PC][2]] = abs(S_noBlocked[Key_FC_OneMet ,FC_OneMet[Key_FC_OneMet][2]])/ abs(S_noBlocked[Key_FC_OneMet ,FC_OneMet[Key_FC_OneMet][1]])

            # V_j = |a/b| * V_i

            Fc_Coefficients[PC[key_PC][2],PC[key_PC][1]] = abs(S_noBlocked[Key_FC_OneMet ,FC_OneMet[Key_FC_OneMet][1]])/ abs(S_noBlocked[Key_FC_OneMet ,FC_OneMet[Key_FC_OneMet][2]])

        else

            # Transpose S:

            S_noBlocked_transpose = S_noBlocked'

            # CO

            CO = zeros(col_noBlocked)

            # Set CO vector:

            CO[PC[key_PC][1]] = -1

            # Removing C:

            S_noBlocked_transpose_removedC = S_noBlocked_transpose[1:end .!= PC[key_PC][2], :]

            rowS, colS = size(S_noBlocked_transpose)
            rowC, colC = size(S_noBlocked_transpose_removedC)

            CO = CO[setdiff(1:end, PC[key_PC][2])]

            # Using Gaussian elimination:

            X = S_noBlocked_transpose_removedC \ CO

            Sol = (S_noBlocked_transpose_removedC*X) - (CO)

            # Determining Fully Coupling:

            S_noBlocked_removedC = S_noBlocked_transpose_removedC'
            if isapprox(norm(Sol), 0.0, atol = epsilon)

                # Editing fctable and changing partially Coupling to Fully Coupling:

                fctable[PC[key_PC][1],PC[key_PC][2]] = 1.0
                fctable[PC[key_PC][2],PC[key_PC][1]] = 1.0

                # Casting SparseVector{Float64, Int64} to SparseMatrixCSC{Float64, Int64}:

                rowC_S_noBlocked_transpose = sparsevec(S_noBlocked_transpose[PC[key_PC][2], :])
                matrix_rowC_S_noBlocked_transpose = sparse(rowC_S_noBlocked_transpose)
                sparseMatrix_rowC = SparseMatrixCSC(matrix_rowC_S_noBlocked_transpose)

                # Calculationg coefficient:
                C = sparseMatrix_rowC' * X

                ## Setting coefficient for i,j:

                # V_i = C * V_j
                Fc_Coefficients[PC[key_PC][1],PC[key_PC][2]] = C[1]

                # V_j = 1/C * V_i
                Fc_Coefficients[PC[key_PC][2],PC[key_PC][1]] = 1 / C[1]

            end
        end
    end

    removeQFCAProcs()

    return fctable, Fc_Coefficients, Dc_Coefficients, n_blocked
end
end
