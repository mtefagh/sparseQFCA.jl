include("../sparseQFCA.jl")
using .sparseQFCA, SparseArrays, DelimitedFiles, Test
S = sparse(readdlm("S.csv", header = false))
@assert typeof(S) == SparseMatrixCSC{Float64,Int64}
rev = readdlm("rev.csv", header = false)[:, 1] .== 1
@assert typeof(rev) == BitArray{1}
certificates, blocked, fctable = @time QFCA(S, rev)
@test all(readdlm("fctable.csv", header = false) .== fctable)