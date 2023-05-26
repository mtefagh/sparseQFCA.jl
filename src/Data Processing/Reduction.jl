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

using .DistributedQFCA, .Pre_processing, COBREXA, SparseArrays

"""
    reduction(myModel)



# INPUTS

- `myModel`:                   A model that has been built using COBREXA's `load_model` function.

# OPTIONAL INPUTS

- `printLevel`:                Verbose level (default: 1). Mute all output with `printLevel = 0`.

# OUTPUTS

-

# EXAMPLES

- Full input/output example
```julia
julia>
```

See also: `dataOfModel()`, `distributedQFCA()`

"""

function reduction(myModel::StandardModel)

    fctable, blocked_index, Fc_Coefficients, Dc_Coefficients = distributedQFCA(myModel)
    row, col = size(fctable)

    A_rows_original = Array{Int64}([])
    A_cols_reduced = Array{Int64}([])
    n = row + length(blocked_index)
    println("Nmber of Original Reactions: $n")
    Reaction_Ids = collect(1:n)
    Reaction_Ids_noBlocked = setdiff(Reaction_Ids, blocked_index)
    A_rows_original = copy(Reaction_Ids)
    A_cols_reduced  = copy(Reaction_Ids)
    println("Reactions Ids : $(length(Reaction_Ids))")
    println(Reaction_Ids)

    ## FC

    FC = Dict()
    FC_Final = Dict()
    FC_Coef = Dict()
    remove_list = []
    c = 1
    for i in range(1, col)
        for j in range(i+1, col)
            if (fctable[i,j] == fctable[j,i] == 1.0) && (i != j)
                FC_Coef[c] = Reaction_Ids_noBlocked[i],Reaction_Ids_noBlocked[j], Fc_Coefficients[i,j]
                FC[c] = Reaction_Ids_noBlocked[i],Reaction_Ids_noBlocked[j]
                c = c + 1
            end
        end
    end

    println("FC:")
    for key in sort(collect(keys(FC)))
        println("$(key) : $(FC[key])")
    end
    println("----------------------------------")
    println("FC with Coef::")
    for key in sort(collect(keys(FC_Coef)))
        println("$(key) : $(FC_Coef[key])")
    end

    s = 1
    for key in sort(collect(keys(FC)))
        if (key >= 2) && (FC[key][1] == FC[key-1][1])
            temp_list = []
            delete_list = []
            for i = 1 : key-1
                if FC[key][1] .== FC[i][1]
                    append!(temp_list, FC[i][2])
                    append!(delete_list, i)
                end
            end
        append!(temp_list, FC[key][2])
        append!(remove_list, temp_list)
        FC_Final[s] = FC[key][1], temp_list
        s = s + 1
        for i in delete_list
            delete!(FC_Final, i)
        end
        else
            FC_Final[s] = FC[s]
            s = s + 1
        end
    end

    remove_list = unique(remove_list)
    remove_list = sort(remove_list)
    println(remove_list)

    for key in sort(collect(keys(FC_Final)))
        if FC_Final[key][1] in remove_list
            delete!(FC_Final, key)
        end
    end

    FC_cluster_members = Array{Int64}([])
    FC_representatives= Array{Int64}([])
    for key in sort(collect(keys(FC_Final)))
        append!(FC_cluster_members, FC_Final[key][2])
        append!(FC_representatives, FC_Final[key][1])
    end

    FC_Clusters = Dict()
    c = 1
    for key in sort(collect(keys(FC_Final)))
        FC_Clusters[c] = FC_Final[key]
        c += 1
    end

    for key in sort(collect(keys(FC_Clusters)))
        println("$(key) : $(FC_Clusters[key])")
    end

    ## DC

    Tolerance = getTolerance()
    remove_list_DC = Array{Int64}([])
    for i in range(1, row)
        DC_row = Dc_Coefficients[i, :]
        DC_row = sparsevec(DC_row)
        if nnz(DC_row) >= 2
            printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)
            println("$i:")
            println(nnz(DC_row))
            for j = 1:length(DC_row)
                if isapprox(DC_row[j] , 0, atol=Tolerance)
                    Dc_Coefficients[i, j] = 0.0
                end
            end
            DC_row = Dc_Coefficients[i, :]
            DC_row = sparsevec(DC_row)
            println(nnz(DC_row))
            for j = 1:length(DC_row)
                if DC_row[j] != 0
                    println("$j : $(DC_row[j])")
                end
            end
            append!(remove_list_DC, i)
        end
    end

    blocked_index = sort(blocked_index)
    FC_cluster_members = sort(FC_cluster_members)
    remove_list_DC = sort(remove_list_DC)

    final = union(blocked_index, FC_cluster_members, remove_list_DC)
    A_cols_reduced = A_cols_reduced[setdiff(range(1, n), final)]

    printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

    println("Blocked : $(length(blocked_index))")
    println(blocked_index)

    printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

    println("FC_Clusters: $(length(FC_cluster_members))")
    for key in sort(collect(keys(FC_Clusters)))
        println("$(key) : $(FC_Clusters[key])")
    end
    println("FC Representatives: $(length(FC_representatives))")
    println(FC_representatives)
    println("FC Cluster Members: $(length(FC_cluster_members))")
    println(FC_cluster_members)

    printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

    println("DCEs LHS: $(length(remove_list_DC))")
    println(remove_list_DC)

    printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

    n = length(A_rows_original)
    n_tylda = length(A_cols_reduced)
    Eliminations = setdiff(A_rows_original, A_cols_reduced)

    println("A Matrix : $n x $n_tylda")
    println("Eliminations: ")
    println("Original: $n")
    println(A_rows_original)
    println("Reduced: $n_tylda")
    println(A_cols_reduced)
    println("Eliminations: $(length(Eliminations))")
    println(Eliminations)

    A = zeros(n, n_tylda)
    println(typeof(A))

    ## Blocked

    for i = 1:n
        if i in blocked_index
            A[i, :] .= 0.0
        end
    end

    ## I

    for i = 1:n_tylda
        A[A_cols_reduced[i],i] = 1.0
    end

    d_1 = sum(A .== 1.0)
    println("Number of 1's in I: $d_1")

    ## FC

    for key in sort(collect(keys(FC_Coef)))

        if FC_Coef[key][1] in A_cols_reduced
            println(FC_Coef[key][1])
            println(FC_Coef[key][2])
            println(FC_Coef[key][3])
            index = findfirst(x -> x == FC_Coef[key][1], A_cols_reduced)
            println(index)
            A[FC_Coef[key][2], index] = FC_Coef[key][3]
            println("-------------------------------------")
        end
    end

    index_1 = findfirst(x -> x == 14, A_cols_reduced)
    index_2 = findfirst(x -> x == 48, A_cols_reduced)

    println("23,$(A_cols_reduced[index_1]) ---> $(A[23,index_1])")
    println("57,$(A_cols_reduced[index_2]) ---> $(A[57,index_2])")
    println("76,$(A_cols_reduced[index_2]) ---> $(A[76,index_2])")

    return A
end

end
