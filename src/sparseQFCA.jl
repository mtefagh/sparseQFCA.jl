#-------------------------------------------------------------------------------------------
#=
    Purpose:    Include all relevant files for running sparseQFCA
    Author:     Iman Ghadimi, Mojtaba Tefagh - Sharif University of Technology
    Date:       July 2023
=#
#-------------------------------------------------------------------------------------------

module sparseQFCA

    using COBREXA, SparseArrays, GLPK, JuMP, LinearAlgebra, Distributed, SharedArrays, SparseArrays

    include("Pre_Processing/Pre_processing.jl")
    include("Pre_Processing/Consistency Checking/TheNaiveApproach.jl")
    include("Pre_Processing/Consistency Checking/SwiftCC.jl")
    include("FCA/distributedQFCA.jl")
    include("FCA/SQFCA.jl")
    include("QuantomeRedNet/Reduction.jl")

    using .Pre_processing, .TheNaiveApproach, .SwiftCC, .DistributedQFCA, .SQFCA, .Reduction

end

#-------------------------------------------------------------------------------------------
