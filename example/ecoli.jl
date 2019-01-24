module ecoli
export S, rev, fctest

using SparseArrays, DelimitedFiles

S = sparse(readdlm("./example/S.csv", header = false))
@assert typeof(S) == SparseMatrixCSC{Float64,Int64}

rev = readdlm("./example/rev.csv", header = false)[:, 1] .== 1
@assert typeof(rev) == BitArray{1}

fctest(fctable) = all(readdlm("./example/fctable.csv", header = false) .== fctable)

end