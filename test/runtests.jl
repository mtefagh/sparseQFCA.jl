
# importing the example model and the sparseQFCA module:

include("ecoli.jl")
using .ecoli, COBREXA, JuMP, sparseQFCA, Test

# finding all the flux coupling relations among the reactions:

fctable = @time QFCA(S, rev)[end]
@test fctest(fctable)

# importing the example model and the TheNaiveApproach library:

include("../example/find_blocked_data.jl")
include("../src/pre_processing.jl")
include("../src/Consistency Checking/TheNaiveApproach.jl")
include("../src/Consistency Checking/SwiftCC.jl")
using .find_blocked_data, .pre_processing, .TheNaiveApproach, .SwiftCC

# Comparing TheNaiveApproach and SwiftCC Outputs:

# iAM_Pv461

blockedList_TheNaive_iAM_Pv461 = @time find_blocked_reactions(myModel_iAM_Pv461)
blockedList_swiftCC_iAM_Pv461 = @time swifCC(myModel_iAM_Pv461)
@test blockedTest_iAM_Pv461(blockedList_TheNaive_iAM_Pv461, blockedList_swiftCC_iAM_Pv461)

# iAT_PLT_636

blockedList_TheNaive_iAT_PLT_636 = @time find_blocked_reactions(myModel_iAT_PLT_636)
blockedList_swiftCC_iAT_PLT_636 = @time swifCC(myModel_iAT_PLT_636)
@test blockedTest_iAT_PLT_636(blockedList_TheNaive_iAT_PLT_636, blockedList_swiftCC_iAT_PLT_636)

# iCN900

blockedList_TheNaive_iCN900 = @time find_blocked_reactions(myModel_iCN900)
blockedList_swiftCC_iCN900 = @time swifCC(myModel_iCN900)
@test blockedTest_iCN900(blockedList_TheNaive_iCN900, blockedList_swiftCC_iCN900)

# iNF517

blockedList_TheNaive_iNF517 = @time find_blocked_reactions(myModel_iNF517)
blockedList_swiftCC_iNF517 = @time swifCC(myModel_iNF517)
@test blockedTest_iNF517(blockedList_TheNaive_iNF517, blockedList_swiftCC_iNF517)

# iNJ661

blockedList_TheNaive_iNJ661 = @time find_blocked_reactions(myModel_iNJ661)
blockedList_swiftCC_iNJ661 = @time swifCC(myModel_iNJ661)
@test blockedTest_iNJ661(blockedList_TheNaive_iNJ661, blockedList_swiftCC_iNJ661)
