module TestData

using Distributed

import JSONFBCModels: JSONFBCModel

import AbstractFBCModels as A

export myModel_e_coli_core, myModel_iIS312, myModel_iAB_RBC_283, fctable_FFCA_e_coli_core, fctable_FFCA_iIS312,
       blockedTest_e_coli_core, blockedTest_iIS312, QFCATest_iIS312, distributedQFCATest_e_coli_core, distributedQFCATest_iIS312

using COBREXA, DelimitedFiles, SparseArrays

## Load Metabolic Networks

myModel_e_coli_core = load_model(JSONFBCModel, "Models/e_coli_core.json", A.CanonicalModel.Model)
myModel_iIS312 = load_model(JSONFBCModel, "Models/iIS312.json", A.CanonicalModel.Model)

## read the CSV file into a matrix of type String

fctable_FFCA_e_coli_core = readdlm("FCTables/MatlabFFCA_Fctable_Ecolicore.csv", ',', Int, header=false)
fctable_FFCA_iIS312 = readdlm("FCTables/MatlabFFCA_Fctable_iIS312.csv", ',', Int, header=false)

## Define functions to Comparing Results

blockedTest_e_coli_core(list_TheNaive, list_SwiftCC) = all(list_TheNaive .== list_SwiftCC)
blockedTest_iIS312(list_TheNaive, list_SwiftCC) = all(list_TheNaive .== list_SwiftCC)
QFCATest_iIS312(fctable_QFCA_iIS312) = all(fctable_QFCA_iIS312 .== fctable_FFCA_iIS312)
distributedQFCATest_e_coli_core(fctable_distributedQFCA_e_coli_core) = all(fctable_distributedQFCA_e_coli_core .== fctable_FFCA_e_coli_core)
distributedQFCATest_iIS312(fctable_distributedQFCA_iIS312) = all(fctable_distributedQFCA_iIS312 .== fctable_FFCA_iIS312)

end
