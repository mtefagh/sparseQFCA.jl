module find_blocked_data

export myModel_ecoli_core, myModel_iAB_RBC_283, myModel_iAF692,# myModel_iCN900, myModel_Recon3D,
       blocked_ecoli_core, blocked_iAB_RBC_283, blocked_iAF692, #blocked_iCN900, blocked_Recon3D,
       blockedTest_ecoli_core, blockedTest_iAB_RBC_283, blockedTest_iAF692,# blockedTest_iCN900, blockedTest_Recon3D
       

using COBREXA, DelimitedFiles

myModel_ecoli_core = load_model(StandardModel,"../example/Models/e_coli_core.xml")
myModel_iAB_RBC_283 = load_model(StandardModel,"../example/Models/iAB_RBC_283.xml")
myModel_iAF692 = load_model(StandardModel,"../example/Models/iAF692.xml")
#myModel_iCN900 = load_model(StandardModel,"../example/Models/iCN900.xml")
#myModel_Recon3D = load_model(StandardModel,"../example/Models/Recon3D.xml")

blocked_ecoli_core = readdlm("../example/PythonResults/e_coli_core.csv", header = false)
blocked_iAB_RBC_283 = readdlm("../example/PythonResults/iAB_RBC_283.csv", header = false)
blocked_iAF692 = readdlm("../example/PythonResults/iAF692.csv", header = false)
#blocked_iCN900 = readdlm("../example/PythonResults/iCN900.csv", header = false)
#try
#       blocked_Recon3D = readdlm("../example/PythonResults/Recon3D.csv", header = false)
#catch e
#       blocked_Recon3D = Array{Any}(undef, 0, 2)
#end

blockedTest_ecoli_core(list) = all(list .== blocked_ecoli_core)
blockedTest_iAB_RBC_283(list) = all(list .== blocked_iAB_RBC_283)
blockedTest_iAF692(list) = all(list .== blocked_iAF692)
#blockedTest_iCN900(list) = all(list .== blocked_iCN900)
#blockedTest_Recon3D(list) = all(list .== blocked_Recon3D)

end
