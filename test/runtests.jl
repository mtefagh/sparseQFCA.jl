include("../src/sparseQFCA.jl")
include("../ecoli.jl")
using .sparseQFCA, .ecoli, Test
fctable = QFCA(S, rev)[end]
@test fctest(fctable)