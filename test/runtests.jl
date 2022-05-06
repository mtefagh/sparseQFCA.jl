# importing the example model and the sparseQFCA module
include("ecoli.jl")
using .ecoli, JuMP, sparseQFCA, Test

# finding all the flux coupling relations among the reactions
fctable = @time QFCA(S, rev)[end]
@test fctest(fctable)

# importing the example model and the TheNaiveApproach library
include("../example/find_blocked_data.jl")
include("../src/TheNaiveApproach.jl")
using .find_blocked_data

# finding all blocked reactions among the reactions

blockedList_e_coli_core = @time find_blocked_reactions(myModel_e_coli_core)[end]
@test blockedTest_e_coli_core(blockedList_e_coli_core)

blockedList_iAB_RBC_283 = @time find_blocked_reactions(myModel_iAB_RBC_283)[end]
@test blockedTest_iAB_RBC_283(blockedList_iAB_RBC_283)

blockedList_iAF692 = @time find_blocked_reactions(myModel_iAF692)[end]
@test blockedTest(blockedList_iAF692)

blockedList_iCN900 = @time find_blocked_reactions(myModel_iCN900)[end]
@test blockedTest(blockedList_iCN900)

blockedList_Recon3D = @time find_blocked_reactions(myModel_Recon3D)[end]
@test blockedTest(blockedList_Recon3D)
