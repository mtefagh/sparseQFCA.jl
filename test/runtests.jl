
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

# 8P

n = 8
addQFCAProcs(n)

# e_coli_core

blockedList_TheNaive_e_coli_core = @time find_blocked_reactions(myModel_e_coli_core)
blockedList_swiftCC_e_coli_core = @time swiftCC(myModel_e_coli_core)
@test blockedTest_e_coli_core(blockedList_TheNaive_e_coli_core, blockedList_swiftCC_e_coli_core)

removeQFCAProcs()

# 8P

n = 8
addQFCAProcs(n)

# iIS312

blockedList_TheNaive_iIS312 = @time find_blocked_reactions(myModel_iIS312)
blockedList_swiftCC_iIS312 = @time swiftCC(myModel_iIS312)
@test blockedTest_iIS312(blockedList_TheNaive_iIS312, blockedList_swiftCC_iIS312)

removeQFCAProcs()

## Comparing final flux coupling table between distributedQFCA and FFCA Algorithms:

# 8P

n = 8
addQFCAProcs(n)

include("TestData.jl")
include("../src/Data Processing/pre_processing.jl")
include("../src/Consistency Checking/SwiftCC.jl")
using .TestData, .pre_processing, .SwiftCC

fctable_distributedQFCA_e_coli_core = @time distributedQFCA(myModel_e_coli_core)
fctable_distributedQFCA_iIS312 = @time distributedQFCA(myModel_iIS312)

removeQFCAProcs()

# Testing

@test distributedQFCATest_e_coli_core(fctable_distributedQFCA_e_coli_core)
@test distributedQFCATest__iIS312(fctable_distributedQFCA_iIS312)
