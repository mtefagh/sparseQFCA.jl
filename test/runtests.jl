### Testing ###

## Adding processes:

using Distributed
addprocs(7)

## Importing Libraries:

include("ecoli.jl")
include("TestData.jl")
include("../src/Data Processing/pre_processing.jl")
include("../src/Consistency Checking/TheNaiveApproach.jl")
@everywhere include("../src/Consistency Checking/SwiftCC.jl")
@everywhere include("../src/QFCA/distributedQFCA.jl")

## sparseQFCA:

using .ecoli, .TestData, .pre_processing, .TheNaiveApproach, .SwiftCC, .DistributedQFCA, COBREXA, JuMP, sparseQFCA, Test, Distributed

# finding all the flux coupling relations among the reactions:

fctable = @time QFCA(S, rev)[end]
@test fctest(fctable)

## Consistency_Checking:

# Comparing the Index of blocked reactions of TheNaiveApproach and SwiftCC:

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

## distributedQFCA:

# Comparing the final flux coupling matrix between distributedQFCA and FFCA Algorithms:

fctable_distributedQFCA_e_coli_core, Fc_Coefficients_e_coli_core, Dc_Coefficients_e_coli_core = @time distributedQFCA(myModel_e_coli_core)
fctable_distributedQFCA_iIS312, Fc_Coefficients_iIS312, Dc_Coefficients_iIS312 = @time distributedQFCA(myModel_iIS312)

removeQFCAProcs()

# Testing

@test distributedQFCATest_e_coli_core(fctable_distributedQFCA_e_coli_core)
@test distributedQFCATest_iIS312(fctable_distributedQFCA_iIS312)
