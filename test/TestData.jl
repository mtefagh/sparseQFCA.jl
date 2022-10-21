
module TestData

export myModel_e_coli_core, myModel_iAM_Pv461, myModel_iAT_PLT_636, myModel_iNJ661,
       distributedQFCATest_e_coli_core, blockedTest_iAM_Pv461, blockedTest_iAT_PLT_636, blockedTest_iNJ661

using COBREXA, DelimitedFiles

# Loading Metabolic Networks:

myModel_e_coli_core = load_model(StandardModel,"Models/e_coli_core.xml")
myModel_iAM_Pv461 = load_model(StandardModel,"Models/iAM_Pv461.xml")
myModel_iAT_PLT_636 = load_model(StandardModel,"Models/iAT_PLT_636.xml")
myModel_iNJ661 = load_model(StandardModel,"Models/iNJ661.xml")

# Defining functions to Compare Results:

blockedTest_iAM_Pv461(list_TheNaive, list_SwiftCC) = all(list_TheNaive .== list_SwiftCC)
blockedTest_iAT_PLT_636(list_TheNaive, list_SwiftCC) = all(list_TheNaive .== list_SwiftCC)
blockedTest_iNJ661(list_TheNaive, list_SwiftCC) = all(list_TheNaive .== list_SwiftCC)
distributedQFCATest_e_coli_core(fctable_seq, fctable_4P, fctable_8P) = all(fctable_seq .== fctable_4P .== fctable_8P)

end
