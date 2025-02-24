cd(@__DIR__)

# Add worker processes to the Julia distributed computing environment:

using Distributed

addprocs(7)
println("Number of Proccess : $(nprocs())")
println("Number of Workers  : $(nworkers())")

### Import Libraries

using COBREXA, JuMP, Test, Distributed, JuMP, HiGHS, Clarabel, JSON, SparseArrays, LinearAlgebra, SharedArrays

# Include the necessary Julia files:
include("TestData.jl")
@everywhere include("../src/sparseQFCA.jl")

# Import required Julia modules:

using .TestData, .sparseQFCA

import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel: Model
import AbstractFBCModels.CanonicalModel: Reaction, Metabolite, Gene, Coupling
import JSONFBCModels: JSONFBCModel

### sparseQFCA:
#=
# Print a message indicating that sparseQFCA is being run on e_coli_core:
printstyled("sparseQFCA :\n"; color=:yellow)
printstyled("iIS312 :\n"; color=:yellow)

# Extracte relevant data from input model:

S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, n_genes_iIS312, lb_iIS312, ub_iIS312, c_vector = sparseQFCA.dataOfModel(myModel_iIS312)
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

S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, m_e_coli_core, n_e_coli_core, n_genes_e_coli_core, lb_e_coli_core, ub_e_coli_core, c_vector = sparseQFCA.dataOfModel(myModel_e_coli_core)
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
S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, n_genes_iIS312, lb_iIS312, ub_iIS312, c_vector = sparseQFCA.dataOfModel(myModel_iIS312)
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
S_e_coli_core, Metabolites_e_coli_core, Reactions_e_coli_core, Genes_e_coli_core, m_e_coli_core, n_e_coli_core, n_genes_e_coli_core, lb_e_coli_core, ub_e_coli_core, c_vector = sparseQFCA.dataOfModel(myModel_e_coli_core)
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
S_iIS312, Metabolites_iIS312, Reactions_iIS312, Genes_iIS312, m_iIS312, n_iIS312, n_genes_iIS312, lb_iIS312, ub_iIS312, c_vector = sparseQFCA.dataOfModel(myModel_iIS312)
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
=#
## ToyModel

ToyModel = Model()
ModelName = "ToyModel"

printstyled("$ModelName :\n"; color=:yellow)

# Genes:
for i = 1:9
    gene = "b" * "$i"
    ToyModel.genes[gene] = Gene()
end

## Metabolites

# IntraCellular:

#m1c
ToyModel.metabolites["m1"] = Metabolite(name = "M1_c", compartment = "inside")
#m2c
ToyModel.metabolites["m2"] = Metabolite(name = "M2_c", compartment = "inside")
#m3c
ToyModel.metabolites["m3"] = Metabolite(name = "M3_c", compartment = "inside")
#m4c
ToyModel.metabolites["m4"] = Metabolite(name = "M4_c", compartment = "inside")

# ExtraCellular:

#m1e
ToyModel.metabolites["m5"] = Metabolite(name = "M1_e", compartment = "outside")
#m3e
ToyModel.metabolites["m6"] = Metabolite(name = "M3_e", compartment = "outside")

## Reactions

M = sparseQFCA.getM(0)

# Forward:

ToyModel.reactions["M1t"] = Reaction(
    name = "transport m1",
    lower_bound = 0.0,
    upper_bound = M,
    stoichiometry = Dict("m5" => -1.0, "m1" => 1.0),
    gene_association_dnf = [["G1","G2"],["G3"]],
    objective_coefficient = 0.0,
)

ToyModel.reactions["rxn2"] = Reaction(
    name = "rxn2",
    lower_bound = 0.0,
    upper_bound = M,
    stoichiometry = Dict("m1" => -2.0, "m2" => 1.0, "m3" => 1.0),
    gene_association_dnf = [["G2"], ["G3"]],
    objective_coefficient = 0.0,
)

ToyModel.reactions["rxn3"] = Reaction(
    name = "rxn3",
    lower_bound = 0.0,
    upper_bound = M,
    stoichiometry = Dict("m2" => -1.0, "m3" => 1.0),
    gene_association_dnf = [["G3","G4"],["G5","G6"]],
    objective_coefficient = 0.0,
)

ToyModel.reactions["M2t"] = Reaction(
    name = "transport m2",
    lower_bound = 0.0,
    upper_bound = M,
    stoichiometry = Dict("m2" => -1.0, "m5" => 1.0),
    gene_association_dnf = [["G4"], ["G1","G7"], ["G3","G5"]],
    objective_coefficient = 0.0,
)

# Foward and Backward:

ToyModel.reactions["rxn1"] = Reaction(
    name = "rxn1",
    lower_bound = -M,
    upper_bound = M,
    stoichiometry = Dict("m1" => -1.0, "m4" => 1.0),
    gene_association_dnf = [["G9"]],
    objective_coefficient = 0.0,
)

ToyModel.reactions["M3t"] = Reaction(
    name = "transport m3",
    lower_bound = -M,
    upper_bound = M,
    stoichiometry = Dict("m3" => -1.0, "m6" => 1.0),
    gene_association_dnf = [["G6"]],
    objective_coefficient = 0.0,
)

# Exchange:

ToyModel.reactions["EX_1"] = Reaction(
    name = "exchange m5",
    lower_bound = -M,
    upper_bound = M,
    stoichiometry = Dict("m5" => -1.0),
    gene_association_dnf = [["G7"]],
    objective_coefficient = 0.0,
)

ToyModel.reactions["EX_2"] = Reaction(
    name = "exchange m6",
    lower_bound = -20,
    upper_bound = M,
    stoichiometry = Dict("m6" => -1.0),
    gene_association_dnf = [["G8"]],
    objective_coefficient = 1.0,
)

ModelName = "ToyModel"  # Define the model name as a string
ToyModel_json = convert(JSONFBCModel, ToyModel)
save_model(ToyModel_json, "../test/Models/$ModelName.json")  # Use the string in the file path
# Read the JSON file
data = JSON.parsefile("Models/$ModelName.json")

# Process reactions to replace '&&' with 'and' and '||' with 'or' in gene_reaction_rule
if haskey(data, "reactions")
    for reaction in data["reactions"]
        if haskey(reaction, "gene_reaction_rule") && !isempty(reaction["gene_reaction_rule"])
            reaction["gene_reaction_rule"] = replace(reaction["gene_reaction_rule"], "&&" => "and", "||" => "or")
        end
    end
end

# Write the corrected JSON file
open("Models/$ModelName.json", "w") do file
    JSON.print(file, data, 1)  # Use 'indent=1' for indentation
end

S_ToyModel, Metabolites_ToyModel, Reactions_ToyModel, Genes_ToyModel, m_ToyModel, n_ToyModel, n_genes_ToyModel, lb_ToyModel, ub_ToyModel, c_vector_ToyModel = sparseQFCA.dataOfModel(ToyModel)

printstyled("FBA - $ModelName:\n"; color=:blue)

# Define the model
FBA_model, solver = sparseQFCA.changeSparseQFCASolver("HiGHS")
# Add decision variables
n = length(Reactions_ToyModel)
@variable(FBA_model, lb_ToyModel[i] <= x[i = 1:n_ToyModel] <= ub_ToyModel[i])
# Set the objective function
@objective(FBA_model, Max, (c_vector_ToyModel)'* x)
@constraint(FBA_model, (S_ToyModel) * x .== 0)
# Solve the model
optimize!(FBA_model)
V_initial = Array{Float64}([])
for i in 1:length(x)
    append!(V_initial, value(x[i]))
end
println(V_initial)

index_c_Toymodel = findfirst(x -> x == 1.0, c_vector_ToyModel)
println("termination_status = $(termination_status(FBA_model))")
println("objective_value = $(objective_value(FBA_model))")
println("Biomass = $(Reactions_ToyModel[index_c_Toymodel]), Flux = $(V_initial[index_c_Toymodel])")

# Ensure that the bounds of all reactions are homogenous
lb_ToyModel, ub_ToyModel = sparseQFCA.homogenization(lb_ToyModel, ub_ToyModel)
# Separate reactions into reversible and irreversible sets:
# Create an array of reaction IDs:
Reaction_Ids_ToyModel = collect(1:n_ToyModel)
irreversible_reactions_id_ToyModel, reversible_reactions_id_ToyModel = sparseQFCA.reversibility(lb_ToyModel, Reaction_Ids_ToyModel)
# Create a new instance of the input model with homogenous bounds:
ModelObject_CC_ToyModel = sparseQFCA.Model_CC(S_ToyModel, Metabolites_ToyModel, Reactions_ToyModel, Genes_ToyModel, m_ToyModel, n_ToyModel, lb_ToyModel, ub_ToyModel)
blocked_index_ToyModel, dualVar_ToyModel = sparseQFCA.swiftCC(ModelObject_CC_ToyModel)
blocked_index_rev_ToyModel = blocked_index_ToyModel ∩ reversible_reactions_id_ToyModel
# Convert to Vector{Int64}:
blocked_index_rev_ToyModel = convert(Vector{Int64}, blocked_index_rev_ToyModel)
# Correct Reversibility:
ModelObject_Crrection_ToyModel = sparseQFCA.Model_Correction(S_ToyModel, Metabolites_ToyModel, Reactions_ToyModel, Genes_ToyModel, m_ToyModel, n_ToyModel, lb_ToyModel, ub_ToyModel, irreversible_reactions_id_ToyModel, reversible_reactions_id_ToyModel)
# Apply distributedReversibility_Correction() to the model and update Reversibility, S and bounds:
SolverName = "HiGHS"
S_ToyModel, lb_ToyModel, ub_ToyModel, irreversible_reactions_id_ToyModel, reversible_reactions_id_ToyModel = sparseQFCA.distributedReversibility_Correction(ModelObject_Crrection_ToyModel, blocked_index_rev_ToyModel, SolverName, false)

printstyled("CorrectedFBA - $ModelName:\n"; color=:blue)

# Define the model
FBA_model_correction, solver = sparseQFCA.changeSparseQFCASolver("HiGHS")
# Add decision variables
@variable(FBA_model_correction, lb_ToyModel[i] <= x[i = 1:n_ToyModel] <= ub_ToyModel[i])
# Set the objective function
println("Objective Function = C'x = $((c_vector_ToyModel)'* x)")
@objective(FBA_model_correction, Max, (c_vector_ToyModel)'* x)
@constraint(FBA_model_correction, (S_ToyModel) * x .== 0)
# Solve the model
optimize!(FBA_model_correction)
V_correction = Array{Float64}([])
for i in 1:length(x)
    append!(V_correction, value(x[i]))
end

println("V_correction:")
println(V_correction)

println("termination_status = $(termination_status(FBA_model_correction))")
println("objective_value = $(objective_value(FBA_model_correction))")
println("Biomass = $(Reactions_ToyModel[index_c_Toymodel]), Flux = $(V_correction[index_c_Toymodel])")

## QuantomeRedNet

ModelName = "ToyModel"
myModel_ToyModel = load_model(JSONFBCModel, "Models/$ModelName.json", A.CanonicalModel.Model)

representatives = []
representatives = convert(Vector{Int64}, representatives)

printstyled("QuantomeRedNet - $ModelName :\n"; color=:yellow)

reducedModelName, A_matrix, reduction_map = sparseQFCA.quantomeReducer(myModel_ToyModel, ModelName, "HiGHS", false, false, representatives)

ToyModel_reduced = load_model(JSONFBCModel, "../src/QuantomeRedNet/ReducedNetworks/$reducedModelName.json", A.CanonicalModel.Model)

S_ToyModelreduced, Metabolites_ToyModelreduced, Reactions_ToyModelreduced, Genes_ToyModelreduced, m_ToyModelreduced, n_ToyModelreduced, n_genes_ToyModelreduced, lb_ToyModelreduced, ub_ToyModelreduced, c_vector_ToyModelreduced = sparseQFCA.dataOfModel(ToyModel_reduced, 0)

Reaction_Ids_ToyModelreduced = collect(1:n_ToyModelreduced)
irreversible_reactions_id_ToyModelreduced, reversible_reactions_id_ToyModelreduced = sparseQFCA.reversibility(lb_ToyModelreduced, Reaction_Ids_ToyModelreduced, 0)

S_ToyModelreduced = dropzeros!(S_ToyModelreduced)

# Find the index of the first occurrence where the element in c_vector is equal to 1.0 in Reduced Network:
index_c_ToyModelreduced = findfirst(x -> x == 1.0, c_vector_ToyModelreduced)

c_vector_ToyModelreduced = A_matrix' * c_vector_ToyModel

S_ToyModelreduced, Metabolites_ToyModelreduced, Metabolites_elimination = sparseQFCA.remove_zeroRows(S_ToyModelreduced, Metabolites_ToyModelreduced)

printstyled("CompressedFBA - $ModelName:\n"; color=:blue)

# Define the model
FBA_model_reduced, solver = sparseQFCA.changeSparseQFCASolver("HiGHS")  # HiGHS supports LP/MILP

# Add decision variables
@variable(FBA_model_reduced, x[1:n_ToyModelreduced])
@variable(FBA_model_reduced, t >= 0)  # Auxiliary variable for norm constraint

t = 1e-3

@constraint(FBA_model_reduced, lb_ToyModel .<= A_matrix * x .<= ub_ToyModel)

# Replace Norm-2 with Norm-Infinity constraints
@constraint(FBA_model_reduced, S_ToyModelreduced * x .<= t)
@constraint(FBA_model_reduced, S_ToyModelreduced * x .>= -t)

# Set objective
@objective(FBA_model_reduced, Max, dot(c_vector_ToyModelreduced, x))

# Solve the model
optimize!(FBA_model_reduced)

# Extract results
V_reduced = []
for i in 1:length(x)
    append!(V_reduced, value(x[i]))
end

println("V_reduced: FBA(Compressed)")
println(V_reduced)

println("termination_status = $(termination_status(FBA_model_reduced))")
println("objective_value = $(objective_value(FBA_model_reduced))")

V = A_matrix * V_reduced

println("V = A * V_reduced")
println(V_initial)

println("Biomass = $(Reactions_ToyModel[index_c_Toymodel]), Flux = $(V_initial[index_c_Toymodel])")

printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:red)

## ToyModel2

printstyled("ToyModel2 :\n"; color=:yellow)

ToyModel2 = Model()

# Genes:
ToyModel2.genes["G1"] = Gene()

## Metabolites

# IntraCellular:

#m1c
ToyModel2.metabolites["m1"] = Metabolite(name = "M1_c", compartment = "inside")

## Reactions

M = sparseQFCA.getM(0)

ToyModel2.reactions["rxn1"] = Reaction(
    name = "rxn1",
    lower_bound = 0.0,
    upper_bound = M,
    stoichiometry = Dict("m1" => 2.0),
    gene_association_dnf = [["G1"]],
    objective_coefficient = 0.0,
)

ToyModel2.reactions["rxn2"] = Reaction(
    name = "rxn2",
    lower_bound = 0.0,
    upper_bound = M,
    stoichiometry = Dict("m1" => -1.0),
    gene_association_dnf = [],
    objective_coefficient = 1.0,
)


ModelName = "ToyModel2"  # Define the model name as a string
ToyModel2_json = convert(JSONFBCModel, ToyModel2)
save_model(ToyModel2_json, "../test/Models/$ModelName.json")  # Use the string in the file path
# Read the JSON file
data = JSON.parsefile("Models/$ModelName.json")

# Process reactions to replace '&&' with 'and' and '||' with 'or' in gene_reaction_rule
if haskey(data, "reactions")
    for reaction in data["reactions"]
        if haskey(reaction, "gene_reaction_rule") && !isempty(reaction["gene_reaction_rule"])
            reaction["gene_reaction_rule"] = replace(reaction["gene_reaction_rule"], "&&" => "and", "||" => "or")
        end
    end
end

# Write the corrected JSON file
open("Models/$ModelName.json", "w") do file
    JSON.print(file, data, 1)  # Use 'indent=1' for indentation
end

S_ToyModel2, Metabolites_ToyModel2, Reactions_ToyModel2, Genes_ToyModel2, m_ToyModel2, n_ToyModel2, n_genes_ToyModel2, lb_ToyModel2, ub_ToyModel2, c_vector_ToyModel2 = sparseQFCA.dataOfModel(ToyModel2)

printstyled("FBA - $ModelName:\n"; color=:blue)

# Define the model
FBA_model, solver = sparseQFCA.changeSparseQFCASolver("HiGHS")
# Add decision variables
n = length(Reactions_ToyModel2)
@variable(FBA_model, lb_ToyModel2[i] <= x[i = 1:n_ToyModel2] <= ub_ToyModel2[i])
# Set the objective function
@objective(FBA_model, Max, (c_vector_ToyModel2)'* x)
@constraint(FBA_model, (S_ToyModel2) * x .== 0)
# Solve the model
optimize!(FBA_model)
V_initial = Array{Float64}([])
for i in 1:length(x)
    append!(V_initial, value(x[i]))
end
println(V_initial)

index_c_ToyModel2 = findfirst(x -> x == 1.0, c_vector_ToyModel2)
println("termination_status = $(termination_status(FBA_model))")
println("objective_value = $(objective_value(FBA_model))")
println("Biomass = $(Reactions_ToyModel2[index_c_ToyModel2]), Flux = $(V_initial[index_c_ToyModel2])")

# Ensure that the bounds of all reactions are homogenous
lb_ToyModel2, ub_ToyModel2 = sparseQFCA.homogenization(lb_ToyModel2, ub_ToyModel2)
# Separate reactions into reversible and irreversible sets:
# Create an array of reaction IDs:
Reaction_Ids_ToyModel2 = collect(1:n_ToyModel2)
irreversible_reactions_id_ToyModel2, reversible_reactions_id_ToyModel2 = sparseQFCA.reversibility(lb_ToyModel2, Reaction_Ids_ToyModel2)
# Create a new instance of the input model with homogenous bounds:
ModelObject_CC_ToyModel2 = sparseQFCA.Model_CC(S_ToyModel2, Metabolites_ToyModel2, Reactions_ToyModel2, Genes_ToyModel2, m_ToyModel2, n_ToyModel2, lb_ToyModel2, ub_ToyModel2)
blocked_index_ToyModel2, dualVar_ToyModel2 = sparseQFCA.swiftCC(ModelObject_CC_ToyModel2)
blocked_index_rev_ToyModel2 = blocked_index_ToyModel2 ∩ reversible_reactions_id_ToyModel2
# Convert to Vector{Int64}:
blocked_index_rev_ToyModel2 = convert(Vector{Int64}, blocked_index_rev_ToyModel2)
# Correct Reversibility:
ModelObject_Crrection_ToyModel2 = sparseQFCA.Model_Correction(S_ToyModel2, Metabolites_ToyModel2, Reactions_ToyModel2, Genes_ToyModel2, m_ToyModel2, n_ToyModel2, lb_ToyModel2, ub_ToyModel2, irreversible_reactions_id_ToyModel2, reversible_reactions_id_ToyModel2)
# Apply distributedReversibility_Correction() to the model and update Reversibility, S and bounds:
SolverName = "HiGHS"
S_ToyModel2, lb_ToyModel2, ub_ToyModel2, irreversible_reactions_id_ToyModel2, reversible_reactions_id_ToyModel2 = sparseQFCA.distributedReversibility_Correction(ModelObject_Crrection_ToyModel2, blocked_index_rev_ToyModel2, SolverName, false)

printstyled("CorrectedFBA - $ModelName:\n"; color=:blue)

# Define the model
FBA_model_correction, solver = sparseQFCA.changeSparseQFCASolver("HiGHS")
# Add decision variables
@variable(FBA_model_correction, lb_ToyModel2[i] <= x[i = 1:n_ToyModel2] <= ub_ToyModel2[i])
# Set the objective function
println("Objective Function = C'x = $((c_vector_ToyModel2)'* x)")
@objective(FBA_model_correction, Max, (c_vector_ToyModel2)'* x)
@constraint(FBA_model_correction, (S_ToyModel2) * x .== 0)
# Solve the model
optimize!(FBA_model_correction)
V_correction = Array{Float64}([])
for i in 1:length(x)
    append!(V_correction, value(x[i]))
end

println("V_correction:")
println(V_correction)

println("termination_status = $(termination_status(FBA_model_correction))")
println("objective_value = $(objective_value(FBA_model_correction))")
println("Biomass = $(Reactions_ToyModel2[index_c_ToyModel2]), Flux = $(V_correction[index_c_ToyModel2])")

## QuantomeRedNet

ModelName = "ToyModel2"
myModel_ToyModel2 = load_model(JSONFBCModel, "Models/$ModelName.json", A.CanonicalModel.Model)

representatives = []
representatives = convert(Vector{Int64}, representatives)

printstyled("QuantomeRedNet - $ModelName :\n"; color=:yellow)

reducedModelName, A_matrix, reduction_map = sparseQFCA.quantomeReducer(myModel_ToyModel2, ModelName, "HiGHS", false, false, representatives)

ToyModel2_reduced = load_model(JSONFBCModel, "../src/QuantomeRedNet/ReducedNetworks/$reducedModelName.json", A.CanonicalModel.Model)

S_ToyModel2reduced, Metabolites_ToyModel2reduced, Reactions_ToyModel2reduced, Genes_ToyModel2reduced, m_ToyModel2reduced, n_ToyModel2reduced, n_genes_ToyModel2reduced, lb_ToyModel2reduced, ub_ToyModel2reduced, c_vector_ToyModel2reduced = sparseQFCA.dataOfModel(ToyModel2_reduced, 0)

Reaction_Ids_ToyModel2reduced = collect(1:n_ToyModel2reduced)
irreversible_reactions_id_ToyModel2reduced, reversible_reactions_id_ToyModel2reduced = sparseQFCA.reversibility(lb_ToyModel2reduced, Reaction_Ids_ToyModel2reduced, 0)

S_ToyModel2reduced = dropzeros!(S_ToyModel2reduced)

# Find the index of the first occurrence where the element in c_vector is equal to 1.0 in Reduced Network:
index_c_ToyModel2reduced = findfirst(x -> x == 1.0, c_vector_ToyModel2reduced)

c_vector_ToyModel2reduced = A_matrix' * c_vector_ToyModel2

S_ToyModel2reduced, Metabolites_ToyModel2reduced, Metabolites_elimination = sparseQFCA.remove_zeroRows(S_ToyModel2reduced, Metabolites_ToyModel2reduced)

printstyled("CompressedFBA - $ModelName:\n"; color=:blue)

# Define the model
FBA_model_reduced, solver = sparseQFCA.changeSparseQFCASolver("HiGHS")  # HiGHS supports LP/MILP

# Add decision variables
@variable(FBA_model_reduced, x[1:n_ToyModel2reduced])
@variable(FBA_model_reduced, t >= 0)  # Auxiliary variable for norm constraint

t = 1e-3

@constraint(FBA_model_reduced, lb_ToyModel2 .<= A_matrix * x .<= ub_ToyModel2)

# Replace Norm-2 with Norm-Infinity constraints
@constraint(FBA_model_reduced, S_ToyModel2reduced * x .<= t)
@constraint(FBA_model_reduced, S_ToyModel2reduced * x .>= -t)

# Set objective
@objective(FBA_model_reduced, Max, dot(c_vector_ToyModel2reduced, x))

# Solve the model
optimize!(FBA_model_reduced)

# Extract results
V_reduced = []
for i in 1:length(x)
    append!(V_reduced, value(x[i]))
end

println("V_reduced: FBA(Compressed)")
println(V_reduced)

println("termination_status = $(termination_status(FBA_model_reduced))")
println("objective_value = $(objective_value(FBA_model_reduced))")

V = A_matrix * V_reduced

println("V = A * V_reduced")
println(V_initial)

println("Biomass = $(Reactions_ToyModel2[index_c_ToyModel2]), Flux = $(V_initial[index_c_ToyModel2])")

printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:red)


#=

# ToyModel3

ToyModel3 = Model()
ModelName = "ToyModel3"

printstyled("$ModelName :\n"; color=:yellow)

# Genes:
for i = 1:5
    gene = "G" * "$i"
    ToyModel3.genes[gene] = Gene()
end

## Metabolites

# IntraCellular:

#m1c
ToyModel3.metabolites["m1"] = Metabolite(name = "M1_c", compartment = "inside")

#m2c
ToyModel3.metabolites["m2"] = Metabolite(name = "M2_c", compartment = "inside")

#m3c
ToyModel3.metabolites["m3"] = Metabolite(name = "M3_c", compartment = "inside")

## Reactions

M = sparseQFCA.getM(0)

# Forward:

ToyModel3.reactions["R1"] = Reaction(
    name = "rxn1",
    lower_bound = -10,
    upper_bound = M,
    stoichiometry = Dict("m1" => -1.0),
    gene_association_dnf = [["G1"]],
    objective_coefficient = 0.0,
)

ToyModel3.reactions["R2"] = Reaction(
    name = "rxn2",
    lower_bound = -M,
    upper_bound = M,
    stoichiometry = Dict("m1" => -2.0, "m2" => 1.0),
    gene_association_dnf = [["G2"]],
    objective_coefficient = 0.0,
)

ToyModel3.reactions["R3"] = Reaction(
    name = "rxn3",
    lower_bound = 0.0,
    upper_bound = M,
    stoichiometry = Dict("m2" => -1.0),
    gene_association_dnf = [["G3"]],
    objective_coefficient = 0.0,
)

ToyModel3.reactions["R4"] = Reaction(
    name = "rxn4",
    lower_bound = 0.0,
    upper_bound = M,
    stoichiometry = Dict("m2" => -1.0, "m3" => 1.0),
    gene_association_dnf = [["G4"]],
    objective_coefficient = 0.0,
)

ToyModel3.reactions["R5"] = Reaction(
    name = "rxn5",
    lower_bound = -M,
    upper_bound = M,
    stoichiometry = Dict("m3" => -1.0),
    gene_association_dnf = [["G5"]],
    objective_coefficient = 1.0,
)

ModelName = "ToyModel3"  # Define the model name as a string
ToyModel3_json = convert(JSONFBCModel, ToyModel3)
save_model(ToyModel3_json, "../test/Models/$ModelName.json")  # Use the string in the file path
# Read the JSON file
data = JSON.parsefile("Models/$ModelName.json")

# Process reactions to replace '&&' with 'and' and '||' with 'or' in gene_reaction_rule
if haskey(data, "reactions")
    for reaction in data["reactions"]
        if haskey(reaction, "gene_reaction_rule") && !isempty(reaction["gene_reaction_rule"])
            reaction["gene_reaction_rule"] = replace(reaction["gene_reaction_rule"], "&&" => "and", "||" => "or")
        end
    end
end

# Write the corrected JSON file
open("Models/$ModelName.json", "w") do file
    JSON.print(file, data, 1)  # Use 'indent=1' for indentation
end

println("FBA - $ModelName:")

S_ToyModel, Metabolites_ToyModel, Reactions_ToyModel, Genes_ToyModel, m_ToyModel, n_ToyModel, n_genes_ToyModel, lb_ToyModel, ub_ToyModel, c_vector_ToyModel = sparseQFCA.dataOfModel(ToyModel3)

# Define the model
FBA_model, solver = sparseQFCA.changeSparseQFCASolver("HiGHS")
# Add decision variables
n = length(Reactions_ToyModel)
@variable(FBA_model, lb_ToyModel[i] <= x[i = 1:n_ToyModel] <= ub_ToyModel[i])
# Set the objective function
@objective(FBA_model, Max, (c_vector_ToyModel)'* x)
@constraint(FBA_model, (S_ToyModel) * x .== 0)
# Solve the model
optimize!(FBA_model)
V = Array{Float64}([])
for i in 1:length(x)
    append!(V, value(x[i]))
end
println("Reactions_Reduced")
println(Reactions_ToyModel)
println(V)

index_c = findfirst(x -> x == 1.0, c_vector_ToyModel)
println("termination_status = $(termination_status(FBA_model))")
println("objective_value = $(objective_value(FBA_model))")
println("Biomass = $(Reactions_ToyModel[index_c]), Flux = $(V[index_c])")

## QuantomeRedNet

ModelName = "ToyModel3"
myModel_ToyModel3 = load_model(JSONFBCModel, "Models/$ModelName.json", A.CanonicalModel.Model)

printstyled("QuantomeRedNet - False - $ModelName :\n"; color=:yellow)
reducedModelName, A_matrix, reduction_map = sparseQFCA.quantomeReducer(myModel_ToyModel3, ModelName, "HiGHS", false, false, representatives)

printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

# ToyModel4

ToyModel4 = Model()
ModelName = "ToyModel4"

printstyled("$ModelName :\n"; color=:yellow)

# Genes:
for i = 1:7
    gene = "G" * "$i"
    ToyModel4.genes[gene] = Gene()
end

## Metabolites

# IntraCellular:

#m1c
ToyModel4.metabolites["m1"] = Metabolite(name = "M1_c", compartment = "inside")

# ExtraCellular:

#m1e
ToyModel4.metabolites["m2"] = Metabolite(name = "M1_e", compartment = "outside")
#m2e
ToyModel4.metabolites["m3"] = Metabolite(name = "M2_e", compartment = "outside")
#m3e
ToyModel4.metabolites["m4"] = Metabolite(name = "M3_e", compartment = "outside")

## Reactions

M = sparseQFCA.getM(0)

ToyModel4.reactions["M1t"] = Reaction(
    name = "transport m1",
    lower_bound = 0.0,
    upper_bound = M,
    stoichiometry = Dict("m2" => -1.0, "m1" => 1.0),
    gene_association_dnf = [["G1"]],
    objective_coefficient = 0.0,
)

ToyModel4.reactions["M2t"] = Reaction(
    name = "transport m2",
    lower_bound = 0.0,
    upper_bound = M,
    stoichiometry = Dict("m1" => -1.0, "m3" => 1.0),
    gene_association_dnf = [["G2"]],
    objective_coefficient = 0.0,
)

ToyModel4.reactions["M3t"] = Reaction(
    name = "transport m3",
    lower_bound = 0.0,
    upper_bound = M,
    stoichiometry = Dict("m1" => -1.0, "m4" => 1.0),
    gene_association_dnf = [["G3"]],
    objective_coefficient = 0.0,
)

# Exchange:

ToyModel4.reactions["EX_1"] = Reaction(
    name = "exchange M1e",
    lower_bound = -5,
    upper_bound = M,
    stoichiometry = Dict("m2" => -1.0),
    gene_association_dnf = [["G4"]],
    objective_coefficient = 0.0,
)

ToyModel4.reactions["EX_2"] = Reaction(
    name = "exchange M2e",
    lower_bound = -7,
    upper_bound = M,
    stoichiometry = Dict("m3" => -1.0),
    gene_association_dnf = [["G5"]],
    objective_coefficient = 0.0,
)

ToyModel4.reactions["EX_3"] = Reaction(
    name = "exchange M3e",
    lower_bound = -M,
    upper_bound = M,
    stoichiometry = Dict("m4" => -1.0),
    gene_association_dnf = [["G6"]],
    objective_coefficient = 1.0,
)

ToyModel4.reactions["rxn1"] = Reaction(
    name = "rxn1",
    lower_bound = 0.0,
    upper_bound = M,
    stoichiometry = Dict("m3" => -1.0, "m2" => 1.0),
    gene_association_dnf = [["G7"]],
    objective_coefficient = 0.0,
)


ModelName = "ToyModel4"  # Define the model name as a string
ToyModel4_json = convert(JSONFBCModel, ToyModel4)
save_model(ToyModel4_json, "../test/Models/$ModelName.json")  # Use the string in the file path
# Read the JSON file
data = JSON.parsefile("Models/$ModelName.json")

# Process reactions to replace '&&' with 'and' and '||' with 'or' in gene_reaction_rule
if haskey(data, "reactions")
    for reaction in data["reactions"]
        if haskey(reaction, "gene_reaction_rule") && !isempty(reaction["gene_reaction_rule"])
            reaction["gene_reaction_rule"] = replace(reaction["gene_reaction_rule"], "&&" => "and", "||" => "or")
        end
    end
end

# Write the corrected JSON file
open("Models/$ModelName.json", "w") do file
    JSON.print(file, data, 1)  # Use 'indent=1' for indentation
end

println("FBA - $ModelName:")

S_ToyModel, Metabolites_ToyModel, Reactions_ToyModel, Genes_ToyModel, m_ToyModel, n_ToyModel, n_genes_ToyModel, lb_ToyModel, ub_ToyModel, c_vector_ToyModel = sparseQFCA.dataOfModel(ToyModel4)

# Define the model
FBA_model, solver = sparseQFCA.changeSparseQFCASolver("HiGHS")
# Add decision variables
n = length(Reactions_ToyModel)
@variable(FBA_model, lb_ToyModel[i] <= x[i = 1:n_ToyModel] <= ub_ToyModel[i])
# Set the objective function
@objective(FBA_model, Max, (c_vector_ToyModel)'* x)
@constraint(FBA_model, (S_ToyModel) * x .== 0)
# Solve the model
optimize!(FBA_model)
V = Array{Float64}([])
for i in 1:length(x)
    append!(V, value(x[i]))
end
println("Reactions_Reduced")
println(Reactions_ToyModel)
println(V)

index_c = findfirst(x -> x == 1.0, c_vector_ToyModel)
println("termination_status = $(termination_status(FBA_model))")
println("objective_value = $(objective_value(FBA_model))")
println("Biomass = $(Reactions_ToyModel[index_c]), Flux = $(V[index_c])")

## QuantomeRedNet

ModelName = "ToyModel4"
myModel_ToyModel4 = load_model(JSONFBCModel, "Models/$ModelName.json", A.CanonicalModel.Model)

printstyled("QuantomeRedNet - False - $ModelName :\n"; color=:yellow)
reducedModelName, A_matrix, reduction_map = sparseQFCA.quantomeReducer(myModel_ToyModel4, ModelName, "HiGHS", false, false, representatives)

printstyled("#-------------------------------------------------------------------------------------------#\n"; color=:yellow)

=#
