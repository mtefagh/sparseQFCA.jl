
# importing the example model and the sparseQFCA module:

include("ecoli.jl")
using .ecoli, COBREXA, JuMP, sparseQFCA, Test

# finding all the flux coupling relations among the reactions:

fctable = @time QFCA(S, rev)[end]
@test fctest(fctable)

# importing the example model and the TheNaiveApproach library:

include("find_blocked_data.jl")
include("../src/pre_processing.jl")
include("../src/Consistency Checking/TheNaiveApproach.jl")
include("../src/Consistency Checking/SwiftCC.jl")
using .find_blocked_data, .pre_processing, .TheNaiveApproach, .SwiftCC

# Comparing TheNaiveApproach and SwiftCC Outputs:

# iAM_Pv461

blockedList_TheNaive_iAM_Pv461 = @time find_blocked_reactions(myModel_iAM_Pv461)
blockedList_swiftCC_iAM_Pv461 = @time swiftCC(myModel_iAM_Pv461)
@test blockedTest_iAM_Pv461(blockedList_TheNaive_iAM_Pv461, blockedList_swiftCC_iAM_Pv461)

# iAT_PLT_636

blockedList_TheNaive_iAT_PLT_636 = @time find_blocked_reactions(myModel_iAT_PLT_636)
blockedList_swiftCC_iAT_PLT_636 = @time swiftCC(myModel_iAT_PLT_636)
@test blockedTest_iAT_PLT_636(blockedList_TheNaive_iAT_PLT_636, blockedList_swiftCC_iAT_PLT_636)

# iNJ661

blockedList_TheNaive_iNJ661 = @time find_blocked_reactions(myModel_iNJ661)
blockedList_swiftCC_iNJ661 = @time swiftCC(myModel_iNJ661)
@test blockedTest_iNJ661(blockedList_TheNaive_iNJ661, blockedList_swiftCC_iNJ661)
