### Test ###

using Distributed

# Add worker processes to the Julia distributed computing environment:
addprocs(7)

## Import Libraries:

include("ecoli.jl")
include("TestData.jl")
include("../src/Data Processing/pre_processing.jl")
include("../src/Consistency Checking/TheNaiveApproach.jl")
@everywhere include("../src/Consistency Checking/SwiftCC.jl")
@everywhere include("../src/QFCA/distributedQFCA.jl")

# Import required Julia modules
using .ecoli, .TestData, .pre_processing, .TheNaiveApproach, .SwiftCC, .DistributedQFCA, COBREXA, JuMP, sparseQFCA, Test, Distributed

## sparseQFCA:

# finding all the flux coupling relations among the reactions:
printstyled("sparseQFCA :\n"; color=:yellow)
printstyled("e_coli_core :\n"; color=:yellow)
fctable = @time QFCA(S, rev)[end]
@test fctest(fctable)
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## Consistency_Checking:

## Compare the Index of blocked reactions of TheNaiveApproach and SwiftCC

## e_coli_core

printstyled("CC_TheNaiveApproach :\n"; color=:yellow)
printstyled("e_coli_core :\n"; color=:yellow)
blockedList_TheNaive_e_coli_core = @time find_blocked_reactions(myModel_e_coli_core)
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)
printstyled("CC_SwiftCC :\n"; color=:yellow)
printstyled("e_coli_core :\n"; color=:yellow)
S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, m_e_coli_core, n_e_coli_core, lb_e_coli_core, ub_e_coli_core = dataOfModel(myModel_e_coli_core)
check_duplicate = check_duplicate_reactions(Reactions_e_coli_core)
lb_e_coli_core, ub_e_coli_core = homogenization(lb_e_coli_core, ub_e_coli_core)
ModelObject_e_coli_core = MyModel(S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, m_e_coli_core, n_e_coli_core, lb_e_coli_core, ub_e_coli_core)
blockedList_swiftCC_e_coli_core, dualVar_e_coli_core = @time swiftCC(ModelObject_e_coli_core)
@test blockedTest_e_coli_core(blockedList_TheNaive_e_coli_core, blockedList_swiftCC_e_coli_core)
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## iIS312

printstyled("CC_TheNaiveApproach :\n"; color=:yellow)
printstyled("iIS312 :\n"; color=:yellow)
blockedList_TheNaive_iIS312 = @time find_blocked_reactions(myModel_iIS312)
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)
printstyled("CC_SwiftCC :\n"; color=:yellow)
printstyled("iIS312 :\n"; color=:yellow)
S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, lb_iIS312, ub_iIS312 = dataOfModel(myModel_iIS312)
check_duplicate = check_duplicate_reactions(Reactions_iIS312)
lb_iIS312, ub_iIS312 = homogenization(lb_iIS312, ub_iIS312)
ModelObject_iIS312 = MyModel(S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, lb_iIS312, ub_iIS312)
blockedList_swiftCC_iIS312, dualVar_e_coli_core_iIS312  = @time swiftCC(ModelObject_iIS312)
@test blockedTest_iIS312(blockedList_TheNaive_iIS312, blockedList_swiftCC_iIS312)
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## distributedQFCA:

## Compare the final flux coupling matrix obtained using two different algorithms, namely distributedQFCA and FFCA

## e_coli_core

printstyled("distributedQFCA :\n"; color=:yellow)
printstyled("e_coli_core :\n"; color=:yellow)
fctable_distributedQFCA_e_coli_core, Fc_Coefficients_e_coli_core, Dc_Coefficients_e_coli_core = @time distributedQFCA(myModel_e_coli_core)
@test distributedQFCATest_e_coli_core(fctable_distributedQFCA_e_coli_core)
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## iIS312

printstyled("distributedQFCA :\n"; color=:yellow)
printstyled("iIS312 :\n"; color=:yellow)
fctable_distributedQFCA_iIS312, Fc_Coefficients_iIS312, Dc_Coefficients_iIS312 = @time distributedQFCA(myModel_iIS312,true)
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)
