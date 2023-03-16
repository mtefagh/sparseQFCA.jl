
module TestData

using Distributed

export myModel_e_coli_core, myModel_iIS312, fctable_FFCA_e_coli_core, fctable_FFCA_iIS312,
       blockedTest_e_coli_core, blockedTest_iIS312, distributedQFCATest_e_coli_core, distributedQFCATest_iIS312

using COBREXA, DelimitedFiles

# Loading Metabolic Networks:

myModel_e_coli_core = load_model(StandardModel,"Models/e_coli_core.xml")
myModel_iIS312 = load_model(StandardModel,"Models/iIS312.xml")
fctable_FFCA_e_coli_core = readdlm("FCTables/MatlabFFCA_Fctable_Ecolicore.csv", header = false)
fctable_FFCA_iIS312 = readdlm("FCTables/MatlabFFCA_Fctable_iIS312.csv", header = false)

# Defining functions to Comparing Results:

blockedTest_e_coli_core(list_TheNaive, list_SwiftCC) = all(list_TheNaive .== list_SwiftCC)
blockedTest_iIS312(list_TheNaive, list_SwiftCC) = all(list_TheNaive .== list_SwiftCC)
distributedQFCATest_e_coli_core(fctable_distributedQFCA_e_coli_core) = all(fctable_distributedQFCA_e_coli_core .== fctable_FFCA_e_coli_core)
distributedQFCATest_iIS312(fctable_distributedQFCA_iIS312) = all(fctable_distributedQFCA_iIS312 .== fctable_FFCA_iIS312)

end
