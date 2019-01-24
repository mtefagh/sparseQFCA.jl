include("../src/sparseQFCA.jl")
include("../example/ecoli.jl")
using .sparseQFCA, .ecoli, Test
fctable = QFCA(S, rev)[end]
@test fctest(fctable)