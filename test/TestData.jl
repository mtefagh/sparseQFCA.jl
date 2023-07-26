module TestData

using Distributed

export myModel_e_coli_core, myModel_iIS312, fctable_FFCA_e_coli_core, fctable_FFCA_iIS312,
       blockedTest_e_coli_core, blockedTest_iIS312, QFCATest_iIS312, distributedQFCATest_e_coli_core, distributedQFCATest_iIS312

using COBREXA, DelimitedFiles, SparseArrays, HTTP


# Create a directory to save the XML models:
output_directory = "Models"
if !isdir(output_directory)
    mkdir(output_directory)
end

## Load Metabolic Networks

## e_coli_core

e_coli_core_model = HTTP.get("http://bigg.ucsd.edu/static/models/e_coli_core.xml")
# Save the e_coli_core_model to a file with the corresponding model name:
model_file_path = joinpath(output_directory, "$e_coli_core.xml")
open(model_file_path, "w") do file
       write(file, e_coli_core_model)
end
myModel_e_coli_core = load_model(StandardModel,"Models/e_coli_core_model.xml")

## iIS312

iIS312_model = HTTP.get("http://bigg.ucsd.edu/static/models/iIS312.xml")
# Save the iIS312_model to a file with the corresponding model name:
model_file_path = joinpath(output_directory, "$iIS312.xml")
open(model_file_path, "w") do file
       write(file, iIS312_model)
end
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
