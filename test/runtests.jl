# importing the example model and the sparseQFCA module
include("ecoli.jl")
using .ecoli, COBREXA, JuMP, sparseQFCA, Test

# finding all the flux coupling relations among the reactions
fctable = @time QFCA(S, rev)[end]
@test fctest(fctable)

# importing the example model and the TheNaiveApproach library

include("../example/find_blocked_data.jl")
include("../src/pre_processing.jl")
include("../src/TheNaiveApproach.jl")

using .find_blocked_data, .pre_processing, .TheNaiveApproach

# finding all blocked reactions among the reactions

blockedList_ecoli_core = @time find_blocked_reactions(myModel_ecoli_core)
@test blockedTest_ecoli_core(blockedList_ecoli_core)

blockedList_iAB_RBC_283 = @time find_blocked_reactions(myModel_iAB_RBC_283)
#@test blockedTest_iAB_RBC_283(blockedList_iAB_RBC_283)

blockedList_iAF692 = @time find_blocked_reactions(myModel_iAF692)
#@test blockedTest_iAF692(blockedList_iAF692)

blockedList_iCN900 = @time find_blocked_reactions(myModel_iCN900)
@test blockedTest_iCN900(blockedList_iCN900)

blockedList_Recon3D = @time find_blocked_reactions(myModel_Recon3D)
@test blockedTest_Recon3D(blockedList_Recon3D)

