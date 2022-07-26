
module find_blocked_data

export myModel_iAM_Pv461, myModel_iAT_PLT_636, myModel_iCN900, myModel_iNF517, myModel_iNJ661,
       blockedTest_iAM_Pv461, blockedTest_iAT_PLT_636, blockedTest_iCN900, blockedTest_iNF517, blockedTest_iNJ661

using COBREXA, DelimitedFiles

# Loading Metabolic Networks:

myModel_iAM_Pv461 = load_model(StandardModel,"../example/Models/iAM_Pv461.xml")
myModel_iAT_PLT_636 = load_model(StandardModel,"../example/Models/iAT_PLT_636.xml")
myModel_iCN900 = load_model(StandardModel,"../example/Models/iCN900.xml")
myModel_iNF517 = load_model(StandardModel,"../example/Models/iNF517.xml")
myModel_iNJ661 = load_model(StandardModel,"../example/Models/iNJ661.xml")

# Defining functions to Compare Results:

blockedTest_iAM_Pv461(list_TheNaive, list_SwiftCC) = all(list_TheNaive .== list_SwiftCC)
blockedTest_iAT_PLT_636(list_TheNaive, list_SwiftCC) = all(list_TheNaive .== list_SwiftCC)
blockedTest_iCN900(list_TheNaive, list_SwiftCC) = all(list_TheNaive .== list_SwiftCC)
blockedTest_iNF517(list_TheNaive, list_SwiftCC) = all(list_TheNaive .== list_SwiftCC)
blockedTest_iNJ661(list_TheNaive, list_SwiftCC) = all(list_TheNaive .== list_SwiftCC)

end
