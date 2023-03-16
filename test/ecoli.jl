module ecoli

export S, rev, fctest, fcfig

using SparseArrays, DelimitedFiles, Colors

S = sparse(readdlm("S.csv", header = false))
@assert typeof(S) == SparseMatrixCSC{Float64,Int64}

rev = readdlm("rev.csv", header = false)[:, 1] .== 1
@assert typeof(rev) == BitArray{1}

fctable = readdlm("fctable.csv", header = false)
fctest(table) = all(table .== fctable)

palette = distinguishable_colors(5)
fcfig = map(x -> palette[convert(Int64, x+1)], fctable)

end
