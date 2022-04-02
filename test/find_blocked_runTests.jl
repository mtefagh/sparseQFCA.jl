
include("../example/find_blocked_data.jl")

using .find_blocked_data, TheNaiveApproach, Test

# finding all blocked reactions among the reactions

blockedList = @time find_blocked_reactions(myModel_e_coli_core)[end]
@test blockedTest(blockedList)

blockedList = @time find_blocked_reactions(myModel_iAB_RBC_283)[end]
@test blockedTest(blockedList)

blockedList = @time find_blocked_reactions(myModel_iAF692)[end]
@test blockedTest(blockedList)

blockedList = @time find_blocked_reactions(myModel_iCN900)[end]
@test blockedTest(blockedList)

blockedList = @time find_blocked_reactions(myModel_Recon3D)[end]
@test blockedTest(blockedList)
