
# importing the example model and the sparseQFCA module:

include("ecoli.jl")
using .ecoli, COBREXA, JuMP, sparseQFCA, Test

# finding all the flux coupling relations among the reactions:

fctable = @time QFCA(S, rev)[end]
@test fctest(fctable)

# importing the example model and the TheNaiveApproach library:

include("TestData.jl")
include("../src/Data Processing/pre_processing.jl")
include("../src/Consistency Checking/TheNaiveApproach.jl")
include("../src/Consistency Checking/SwiftCC.jl")
using .TestData, .pre_processing, .TheNaiveApproach, .SwiftCC

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

## Comparing fctable among 1P and 4P and 8P:

# 1P

include("TestData.jl")
include("../src/Data Processing/pre_processing.jl")
include("../src/Consistency Checking/SwiftCC.jl")
using .TestData, .pre_processing, .SwiftCC

fctable_seq_e_coli_core = @time distributedQFCA(myModel_e_coli_core)

# 4P

n = 4
addQFCAProcs(n)

include("TestData.jl")
include("../src/Data Processing/pre_processing.jl")
include("../src/Consistency Checking/SwiftCC.jl")
using .TestData, .pre_processing, .SwiftCC

fctable_4P_e_coli_core = @time distributedQFCA(myModel_e_coli_core)

removeQFCAProcs()

# 8P

n = 8
addQFCAProcs(n)

include("TestData.jl")
include("../src/Data Processing/pre_processing.jl")
include("../src/Consistency Checking/SwiftCC.jl")
using .TestData, .pre_processing, .SwiftCC

fctable_8P_e_coli_core = @time distributedQFCA(myModel_e_coli_core)

removeQFCAProcs()

# Test

@test distributedQFCATest_e_coli_core(fctable_seq_e_coli_core, fctable_4P_e_coli_core, fctable_8P_e_coli_core)
