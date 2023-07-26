module TestData

using Distributed

export myModel_e_coli_core, myModel_iIS312, fctable_FFCA_e_coli_core, fctable_FFCA_iIS312,
       blockedTest_e_coli_core, blockedTest_iIS312, QFCATest_iIS312, distributedQFCATest_e_coli_core, distributedQFCATest_iIS312

using COBREXA, DelimitedFiles, SparseArrays, HTTP

## Load Metabolic Networks

# e_coli_core:
e_coli_core_model = HTTP.get("http://bigg.ucsd.edu/static/models/e_coli_core.xml")
write("Models/e_coli_core_model.xml",e_coli_core_model.body)
myModel_e_coli_core = load_model(StandardModel,"Models/e_coli_core_model.xml")

# iIS312:
iIS312_model = HTTP.get("http://bigg.ucsd.edu/static/models/iIS312.xml")
write("Models/iIS312_model.xml",ecoli_model.body)
myModel_iIS312 = load_model(StandardModel,"Models/iIS312_model.xml")

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
