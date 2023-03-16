
## sparseQFCA:

include("ecoli.jl")

using .ecoli, COBREXA, JuMP, sparseQFCA, Test, Distributed

# finding all the flux coupling relations among the reactions:

fctable = @time QFCA(S, rev)[end]
@test fctest(fctable)

## Consistency_Checking:

addprocs(7)

include("TestData.jl")
include("../src/Data Processing/pre_processing.jl")
include("../src/Consistency Checking/TheNaiveApproach.jl")
@everywhere include("../src/Consistency Checking/SwiftCC.jl")
@everywhere include("../src/QFCA/distributedQFCA.jl")
using .TestData, .pre_processing, .TheNaiveApproach, .SwiftCC, .DistributedQFCA

# Comparing TheNaiveApproach and SwiftCC Outputs:

# e_coli_core

blockedList_TheNaive_e_coli_core = @time find_blocked_reactions(myModel_e_coli_core)

S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, m_e_coli_core, n_e_coli_core, lb_e_coli_core, ub_e_coli_core = dataOfModel(myModel_e_coli_core)
lb_e_coli_core, ub_e_coli_core = homogenization(lb_e_coli_core, ub_e_coli_core)
ModelObject_e_coli_core = MyModel(S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, m_e_coli_core, n_e_coli_core, lb_e_coli_core, ub_e_coli_core)
blockedList_swiftCC_e_coli_core, dualVar_e_coli_core = @time swiftCC(ModelObject_e_coli_core)
@test blockedTest_e_coli_core(blockedList_TheNaive_e_coli_core, blockedList_swiftCC_e_coli_core)

# iIS312

blockedList_TheNaive_iIS312 = @time find_blocked_reactions(myModel_iIS312)
S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, lb_iIS312, ub_iIS312 = dataOfModel(myModel_iIS312)
lb_iIS312, ub_iIS312 = homogenization(lb_iIS312, ub_iIS312)
ModelObject_iIS312 = MyModel(S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, lb_iIS312, ub_iIS312)
blockedList_swiftCC_iIS312, dualVar_e_coli_core_iIS312  = @time swiftCC(ModelObject_iIS312)
@test blockedTest_iIS312(blockedList_TheNaive_iIS312, blockedList_swiftCC_iIS312)

## QFCA:

# Comparing final flux coupling table between distributedQFCA and FFCA Algorithms:


fctable_distributedQFCA_e_coli_core, Fc_Coefficients_e_coli_core, Dc_Coefficients_e_coli_core = @time distributedQFCA(myModel_e_coli_core)

d_0 = 0
d_1 = 0
d_2 = 0
d_3 = 0
d_4 = 0

d_0 = sum(fctable_distributedQFCA_e_coli_core .== 0.0)
d_1 = sum(fctable_distributedQFCA_e_coli_core .== 1.0)
d_2 = sum(fctable_distributedQFCA_e_coli_core .== 2.0)
d_3 = sum(fctable_distributedQFCA_e_coli_core .== 3.0)
d_4 = sum(fctable_distributedQFCA_e_coli_core .== 4.0)

println()

println("Final fctable_distributedQFCA_e_coli_core : ")

println()

println("Number of 0's (unCoupled) : $d_0")
println("Number of 1's (Fully)     : $d_1")
println("Number of 2's (Partialy)  : $d_2")
println("Number of 3's (DC i-->j)  : $d_3")
println("Number of 4's (DC j-->i)  : $d_4")

fctable_distributedQFCA_e_coli_core, Fc_Coefficients_e_coli_core, Dc_Coefficients_e_coli_core = @time distributedQFCA(myModel_e_coli_core, true)

d_0 = 0
d_1 = 0
d_2 = 0
d_3 = 0
d_4 = 0

d_0 = sum(fctable_distributedQFCA_e_coli_core .== 0.0)
d_1 = sum(fctable_distributedQFCA_e_coli_core .== 1.0)
d_2 = sum(fctable_distributedQFCA_e_coli_core .== 2.0)
d_3 = sum(fctable_distributedQFCA_e_coli_core .== 3.0)
d_4 = sum(fctable_distributedQFCA_e_coli_core .== 4.0)

println()

println("Final fctable_distributedQFCA_e_coli_core : ")

println()

println("Number of 0's (unCoupled) : $d_0")
println("Number of 1's (Fully)     : $d_1")
println("Number of 2's (Partialy)  : $d_2")
println("Number of 3's (DC i-->j)  : $d_3")
println("Number of 4's (DC j-->i)  : $d_4")

fctable_distributedQFCA_iIS312, Fc_Coefficients_iIS312, Dc_Coefficients_iIS312 = @time distributedQFCA(myModel_iIS312)

removeQFCAProcs()

# Testing

@test distributedQFCATest_e_coli_core(fctable_distributedQFCA_e_coli_core)
@test distributedQFCATest__iIS312(fctable_distributedQFCA_iIS312)
