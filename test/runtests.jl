### Test ###

using Distributed

# Add worker processes to the Julia distributed computing environment:
addprocs(7)
println("Number of Proccess : $(nprocs())")
println("Number of Workers  : $(nworkers())")

### Import Libraries

# Include the necessary Julia files:
include("TestData.jl")
@everywhere include("../src/sparseQFCA.jl")

# Import required Julia modules:
using COBREXA, JuMP, Test, Distributed
using .TestData, .sparseQFCA

### sparseQFCA:

# Print a message indicating that sparseQFCA is being run on e_coli_core:
printstyled("sparseQFCA :\n"; color=:yellow)
printstyled("iIS312 :\n"; color=:yellow)

# Extracte relevant data from input model:

S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, Genes_Reactions_iIS312, m_iIS312, n_iIS312, n_genes_iIS312, lb_iIS312, ub_iIS312 = sparseQFCA.dataOfModel(myModel_iIS312)
# Ensure that the bounds of all reactions are homogenous:
lb_iIS312, ub_iIS312 = sparseQFCA.homogenization(lb_iIS312, ub_iIS312)
# Separate reactions into reversible and irreversible sets:
# Create an array of reaction IDs:
Reaction_Ids_iIS312 = collect(1:n_iIS312)
irreversible_reactions_id_iIS312, reversible_reactions_id_iIS312 = sparseQFCA.reversibility(lb_iIS312, Reaction_Ids_iIS312)
# Create a new instance of the input model with homogenous bounds:
ModelObject_CC_iIS312 = sparseQFCA.Model_CC(S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, lb_iIS312, ub_iIS312)
blocked_index_iIS312, dualVar_iIS312 = sparseQFCA.swiftCC(ModelObject_CC_iIS312)
blocked_index_rev_iIS312 = blocked_index_iIS312 ∩ reversible_reactions_id_iIS312
# Convert to Vector{Int64}
blocked_index_rev_iIS312 = convert(Vector{Int64}, blocked_index_rev_iIS312)
# Correct Reversibility:
ModelObject_Crrection_iIS312 = sparseQFCA.Model_Correction(S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, lb_iIS312, ub_iIS312, irreversible_reactions_id_iIS312, reversible_reactions_id_iIS312)
S_iIS312, lb_iIS312, ub_iIS312, irreversible_reactions_id_iIS312, reversible_reactions_id_iIS312 = sparseQFCA.distributedReversibility_Correction(ModelObject_Crrection_iIS312, blocked_index_rev_iIS312)
# Create Rev Vector:
rev_iIS312 = zeros(Bool,n_iIS312)
for i in reversible_reactions_id_iIS312
    rev_iIS312[i] = true
end
# Run QFCA on S and rev, and save the output to fctable:
fctable_QFCA_iIS312 = @time sparseQFCA.QFCA(S_iIS312, rev_iIS312)[end]
# Test that the results of QFCA are correct for the iIS312 model:
@test QFCATest_iIS312(fctable_QFCA_iIS312)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

### Consistency_Checking:

## Compare the Index of blocked reactions of TheNaiveApproach and SwiftCC

## e_coli_core

## Print a message indicating that TheNaiveApproach is being run on e_coli_core

printstyled("CC_TheNaiveApproach :\n"; color=:yellow)
printstyled("e_coli_core :\n"; color=:yellow)
# Find blocked reactions in myModel_e_coli_core using TheNaiveApproach, and save the output to blockedList_TheNaive_e_coli_core:
blockedList_TheNaive_e_coli_core = @time sparseQFCA.find_blocked_reactions(myModel_e_coli_core)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## Print a message indicating that SwiftCC is being run on e_coli_core

printstyled("CC_SwiftCC :\n"; color=:yellow)
printstyled("e_coli_core :\n"; color=:yellow)
# Get the necessary data from myModel_e_coli_core:

S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, Genes_Reactions_e_coli_core, m_e_coli_core, n_e_coli_core, n_genes_e_coli_core, lb_e_coli_core, ub_e_coli_core = sparseQFCA.dataOfModel(myModel_e_coli_core)
# Check for duplicate reactions in Reactions_e_coli_core:
check_duplicate = sparseQFCA.check_duplicate_reactions(Reactions_e_coli_core)
# Homogenize the lower and upper bounds of the reactions in myModel_e_coli_core:
lb_e_coli_core, ub_e_coli_core = sparseQFCA.homogenization(lb_e_coli_core, ub_e_coli_core)
# Create a ModelObject from the data in myModel_e_coli_core:
ModelObject_e_coli_core = sparseQFCA.Model_CC(S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, m_e_coli_core, n_e_coli_core, lb_e_coli_core, ub_e_coli_core)
# Find blocked reactions in the e_coli_core model using the swiftCC method and time the operation:
blockedList_swiftCC_e_coli_core, dualVar_e_coli_core = @time sparseQFCA.swiftCC(ModelObject_e_coli_core)
# Test that the results of the naive approach and swiftCC approach are the same:
@test blockedTest_e_coli_core(blockedList_TheNaive_e_coli_core, blockedList_swiftCC_e_coli_core)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## iIS312

## Print a message indicating that TheNaiveApproach is being run on iIS312

printstyled("CC_TheNaiveApproach :\n"; color=:yellow)
printstyled("iIS312 :\n"; color=:yellow)
# Find blocked reactions in the iIS312 model and time the operation:
blockedList_TheNaive_iIS312 = @time sparseQFCA.find_blocked_reactions(myModel_iIS312)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## Print a message indicating that SwiftCC is being run on iIS312

printstyled("CC_SwiftCC :\n"; color=:yellow)
printstyled("iIS312 :\n"; color=:yellow)

# Get data from the iIS312 model:
S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, Genes_Reactions_iIS312, m_iIS312, n_iIS312, n_genes_iIS312, lb_iIS312, ub_iIS312 = sparseQFCA.dataOfModel(myModel_iIS312)
# Check for duplicate reactions in the iIS312 model:
check_duplicate = sparseQFCA.check_duplicate_reactions(Reactions_iIS312)
# Homogenize the lower and upper bounds for reactions in the iIS312 model:
lb_iIS312, ub_iIS312 = sparseQFCA.homogenization(lb_iIS312, ub_iIS312)
# Create a model object from the iIS312 model data:
ModelObject_iIS312 = sparseQFCA.Model_CC(S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, lb_iIS312, ub_iIS312)
# Find blocked reactions in the iIS312 model using the swiftCC method and time the operation:
blockedList_swiftCC_iIS312, dualVar_e_coli_core_iIS312  = @time sparseQFCA.swiftCC(ModelObject_iIS312)
# Test that the results of the naive approach and swiftCC approach are the same:
@test blockedTest_iIS312(blockedList_TheNaive_iIS312, blockedList_swiftCC_iIS312)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

### distributedQFCA:

## Print a message indicating that distributedQFCA is being run on e_coli_core

printstyled("distributedQFCA :\n"; color=:yellow)
printstyled("e_coli_core :\n"; color=:yellow)

# Extracte relevant data from input model:
S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, Genes_Reactions_e_coli_core, m_e_coli_core, n_e_coli_core, n_genes_e_coli_core, lb_e_coli_core, ub_e_coli_core = sparseQFCA.dataOfModel(myModel_e_coli_core)
# Ensure that the bounds of all reactions are homogenous
lb_e_coli_core, ub_e_coli_core = sparseQFCA.homogenization(lb_e_coli_core, ub_e_coli_core)
# Separate reactions into reversible and irreversible sets:
# Create an array of reaction IDs:
Reaction_Ids_e_coli_core = collect(1:n_e_coli_core)
irreversible_reactions_id_e_coli_core, reversible_reactions_id_e_coli_core = sparseQFCA.reversibility(lb_e_coli_core, Reaction_Ids_e_coli_core)
# Create a new instance of the input model with homogenous bounds:
ModelObject_CC_e_coli_core = sparseQFCA.Model_CC(S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, m_e_coli_core, n_e_coli_core, lb_e_coli_core, ub_e_coli_core)
blocked_index_e_coli_core, dualVar_e_coli_core = sparseQFCA.swiftCC(ModelObject_CC_e_coli_core)
blocked_index_rev_e_coli_core = blocked_index_e_coli_core ∩ reversible_reactions_id_e_coli_core
# Convert to Vector{Int64}:
blocked_index_rev_e_coli_core = convert(Vector{Int64}, blocked_index_rev_e_coli_core)
# Correct Reversibility:
ModelObject_Crrection_e_coli_core = sparseQFCA.Model_Correction(S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, m_e_coli_core, n_e_coli_core, lb_e_coli_core, ub_e_coli_core, irreversible_reactions_id_e_coli_core, reversible_reactions_id_e_coli_core)
S_e_coli_core, lb_e_coli_core, ub_e_coli_core, irreversible_reactions_id_e_coli_core, reversible_reactions_id_e_coli_core = sparseQFCA.distributedReversibility_Correction(ModelObject_Crrection_e_coli_core, blocked_index_rev_e_coli_core)
# Run distributedQFCA method on the model and time the operation:
row_e_coli_core, col_e_coli_core = size(S_e_coli_core)
ModelObject_QFCA_e_coli_core = sparseQFCA.Model_QFCA(S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, row_e_coli_core, col_e_coli_core, lb_e_coli_core, ub_e_coli_core, irreversible_reactions_id_e_coli_core, reversible_reactions_id_e_coli_core)
# Run distributedQFCA method on the e_coli_core model and time the operation:
# Convert to Vector{Int64}:
blocked_index_e_coli_core = convert(Vector{Int64}, blocked_index_e_coli_core)
fctable_distributedQFCA_e_coli_core, Fc_Coefficients_e_coli_core, Dc_Coefficients_e_coli_core = @time sparseQFCA.distributedQFCA(ModelObject_QFCA_e_coli_core, blocked_index_e_coli_core)
# convert the shared matrix to a regular matrix:
fctable_distributedQFCA_e_coli_core = convert(Matrix{Int}, fctable_distributedQFCA_e_coli_core)
# Test that the results of distributedQFCA are correct for the e_coli_core model:
@test distributedQFCATest_e_coli_core(fctable_distributedQFCA_e_coli_core)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## Print a message indicating that distributedQFCA is being run on iIS312

printstyled("distributedQFCA :\n"; color=:yellow)
printstyled("iIS312 :\n"; color=:yellow)

# Extracte relevant data from input model:
S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, Genes_Reactions_iIS312, m_iIS312, n_iIS312, n_genes_iIS312, lb_iIS312, ub_iIS312 = sparseQFCA.dataOfModel(myModel_iIS312)
# Ensure that the bounds of all reactions are homogenous
lb_iIS312, ub_iIS312 = sparseQFCA.homogenization(lb_iIS312, ub_iIS312)
# Separate reactions into reversible and irreversible sets:
# Create an array of reaction IDs:
Reaction_Ids_iIS312 = collect(1:n_iIS312)
irreversible_reactions_id_iIS312, reversible_reactions_id_iIS312 = sparseQFCA.reversibility(lb_iIS312, Reaction_Ids_iIS312)
# Create a new instance of the input model with homogenous bounds:
ModelObject_CC_iIS312 = sparseQFCA.Model_CC(S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, lb_iIS312, ub_iIS312)
blocked_index_iIS312, dualVar_iIS312 = sparseQFCA.swiftCC(ModelObject_CC_iIS312)
blocked_index_rev_iIS312 = blocked_index_iIS312 ∩ reversible_reactions_id_iIS312
# Convert to Vector{Int64}:
blocked_index_rev_iIS312 = convert(Vector{Int64}, blocked_index_rev_iIS312)
# Correct Reversibility:
ModelObject_Crrection_iIS312 = sparseQFCA.Model_Correction(S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, lb_iIS312, ub_iIS312, irreversible_reactions_id_iIS312, reversible_reactions_id_iIS312)
S_iIS312, lb_iIS312, ub_iIS312, irreversible_reactions_id_iIS312, reversible_reactions_id_iIS312 = sparseQFCA.distributedReversibility_Correction(ModelObject_Crrection_iIS312, blocked_index_rev_iIS312)
# Run distributedQFCA method on the model and time the operation:
row_iIS312, col_iIS312 = size(S_iIS312)
ModelObject_QFCA_iIS312 = sparseQFCA.Model_QFCA(S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, row_iIS312, col_iIS312, lb_iIS312, ub_iIS312, irreversible_reactions_id_iIS312, reversible_reactions_id_iIS312)
# Run distributedQFCA method on the iIS312 model and time the operation:
# Convert to Vector{Int64}:
blocked_index_iIS312 = convert(Vector{Int64}, blocked_index_iIS312)
fctable_distributedQFCA_iIS312, Fc_Coefficients_iIS312, Dc_Coefficients_iIS312 = @time sparseQFCA.distributedQFCA(ModelObject_QFCA_iIS312, blocked_index_iIS312)
# convert the shared matrix to a regular matrix:
fctable_distributedQFCA_iIS312 = convert(Matrix{Int}, fctable_distributedQFCA_iIS312)
# Test that the results of distributedQFCA are correct for the iIS312 model:
@test distributedQFCATest_iIS312(fctable_distributedQFCA_iIS312)
# Print a separator:
printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

## QuantomeRedNet

printstyled("QuantomeRedNet :\n"; color=:yellow)
printstyled("e_coli_core :\n"; color=:yellow)
model, A, reduct_map = @time sparseQFCA.quantomeReducer(myModel_e_coli_core)

printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)
