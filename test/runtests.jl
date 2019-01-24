include("sparseQFCA.jl")
using .sparseQFCA, SparseArrays, DelimitedFiles
S = sparse(readdlm("S.csv", header = false))
@assert typeof(S) == SparseMatrixCSC{Float64,Int64}
rev = readdlm("rev.csv", header = false)[:, 1] .== 1
@assert typeof(rev) == BitArray{1}
certificates, blocked, fctable = @time QFCA(S, rev)
println("The answer is $(all(readdlm("fctable.csv", header = false) .== fctable) ? "correct" : "wrong").")