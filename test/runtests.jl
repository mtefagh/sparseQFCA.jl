### Test ###

using Distributed

# Add worker processes to the Julia distributed computing environment:
addprocs(7)

## Import Libraries

# Include the necessary Julia files:
include("ecoli.jl")
include("TestData.jl")
include("../src/Data Processing/pre_processing.jl")
include("../src/Consistency Checking/TheNaiveApproach.jl")
@everywhere include("../src/Consistency Checking/SwiftCC.jl")
@everywhere include("../src/QFCA/distributedQFCA.jl")

# Import required Julia modules:
using .ecoli, .TestData, .pre_processing, .TheNaiveApproach, .SwiftCC, .DistributedQFCA, COBREXA, JuMP, sparseQFCA, Test, Distributed

## sparseQFCA:

# Print a message indicating that sparseQFCA is being run on e_coli_core:
printstyled("sparseQFCA :\n"; color=:yellow)
printstyled("e_coli_core :\n"; color=:yellow)

# Run QFCA on S and rev, and save the output to fctable:
fctable = @time QFCA(S, rev)[end]

# Test the output of QFCA using the fctest function:
@test fctest(fctable)

# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## Consistency_Checking:

## Compare the Index of blocked reactions of TheNaiveApproach and SwiftCC

## e_coli_core

## Print a message indicating that TheNaiveApproach is being run on e_coli_core

printstyled("CC_TheNaiveApproach :\n"; color=:yellow)
printstyled("e_coli_core :\n"; color=:yellow)
# Find blocked reactions in myModel_e_coli_core using TheNaiveApproach, and save the output to blockedList_TheNaive_e_coli_core:
blockedList_TheNaive_e_coli_core = @time find_blocked_reactions(myModel_e_coli_core)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## Print a message indicating that SwiftCC is being run on e_coli_core

printstyled("CC_SwiftCC :\n"; color=:yellow)
printstyled("e_coli_core :\n"; color=:yellow)
# Get the necessary data from myModel_e_coli_core:
S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, m_e_coli_core, n_e_coli_core, lb_e_coli_core, ub_e_coli_core = dataOfModel(myModel_e_coli_core)
# Check for duplicate reactions in Reactions_e_coli_core:
check_duplicate = check_duplicate_reactions(Reactions_e_coli_core)
# Homogenize the lower and upper bounds of the reactions in myModel_e_coli_core:
lb_e_coli_core, ub_e_coli_core = homogenization(lb_e_coli_core, ub_e_coli_core)
# Create a ModelObject from the data in myModel_e_coli_core:
ModelObject_e_coli_core = MyModel(S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, m_e_coli_core, n_e_coli_core, lb_e_coli_core, ub_e_coli_core)
# Find blocked reactions in the e_coli_core model using the swiftCC method and time the operation:
blockedList_swiftCC_e_coli_core, dualVar_e_coli_core = @time swiftCC(ModelObject_e_coli_core)
# Test that the results of the naive approach and swiftCC approach are the same:
@test blockedTest_e_coli_core(blockedList_TheNaive_e_coli_core, blockedList_swiftCC_e_coli_core)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## iIS312

## Print a message indicating that TheNaiveApproach is being run on iIS312

printstyled("CC_TheNaiveApproach :\n"; color=:yellow)
printstyled("iIS312 :\n"; color=:yellow)
# Find blocked reactions in the iIS312 model and time the operation:
blockedList_TheNaive_iIS312 = @time find_blocked_reactions(myModel_iIS312)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## Print a message indicating that SwiftCC is being run on iIS312

printstyled("CC_SwiftCC :\n"; color=:yellow)
printstyled("iIS312 :\n"; color=:yellow)

# Get data from the iIS312 model:
S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, lb_iIS312, ub_iIS312 = dataOfModel(myModel_iIS312)
# Check for duplicate reactions in the iIS312 model:
check_duplicate = check_duplicate_reactions(Reactions_iIS312)
# Homogenize the lower and upper bounds for reactions in the iIS312 model:
lb_iIS312, ub_iIS312 = homogenization(lb_iIS312, ub_iIS312)
# Create a model object from the iIS312 model data:
ModelObject_iIS312 = MyModel(S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, lb_iIS312, ub_iIS312)
# Find blocked reactions in the iIS312 model using the swiftCC method and time the operation:
blockedList_swiftCC_iIS312, dualVar_e_coli_core_iIS312  = @time swiftCC(ModelObject_iIS312)
# Test that the results of the naive approach and swiftCC approach are the same:
@test blockedTest_iIS312(blockedList_TheNaive_iIS312, blockedList_swiftCC_iIS312)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## distributedQFCA:

## Print a message indicating that distributedQFCA is being run on e_coli_core

printstyled("distributedQFCA :\n"; color=:yellow)
printstyled("e_coli_core :\n"; color=:yellow)

# Run distributedQFCA method on the e_coli_core model and time the operation:
fctable_distributedQFCA_e_coli_core, Fc_Coefficients_e_coli_core, Dc_Coefficients_e_coli_core = @time distributedQFCA(myModel_e_coli_core)
# convert the shared matrix to a regular matrix:
fctable_distributedQFCA_e_coli_core = convert(Matrix{Int}, fctable_distributedQFCA_e_coli_core)
# Test that the results of distributedQFCA are correct for the e_coli_core model:
@test distributedQFCATest_e_coli_core(fctable_distributedQFCA_e_coli_core)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## Print a message indicating that distributedQFCA is being run on iIS312

printstyled("distributedQFCA :\n"; color=:yellow)
printstyled("iIS312 :\n"; color=:yellow)
# Run distributedQFCA method on the iIS312 model and time the operation:
fctable_distributedQFCA_iIS312, Fc_Coefficients_iIS312, Dc_Coefficients_iIS312 = @time distributedQFCA(myModel_iIS312,true)
# convert the shared matrix to a regular matrix:
fctable_distributedQFCA_e_coli_core = convert(Matrix{Int}, fctable_distributedQFCA_iIS312)  
# Test that the results of distributedQFCA are correct for the iIS312 model:
@test distributedQFCATest_iIS312(fctable_distributedQFCA_iIS312)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)
