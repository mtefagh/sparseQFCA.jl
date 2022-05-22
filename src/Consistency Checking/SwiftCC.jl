# SwiftCC

using COBREXA
using JuMP
using GLPK
using LinearAlgebra

myModel = load_model(StandardModel,"/Users/iman/Desktop/M.S.Thesis/SourceCode/Data/e_coli_core.xml")

function dataOfModel(myModel)
   S = stoichiometry(myModel)
   Metabolites = metabolites(myModel)
   Reactions = reactions(myModel)
   Genes = genes(myModel)
   m = length(metabolites(myModel))
   n = length(reactions(myModel))
   lb = lower_bounds(myModel)
   ub = upper_bounds(myModel)
   return S, Metabolites, Reactions, Genes, m, n, lb, ub
end

S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel)

function reversibility(myModel)
   n = length(reactions(myModel))
   irreversible_reactions_id = []
   reversible_reactions_id = []
   for i in 1:n
       if lower_bounds(myModel)[i] >= 0
           append!(irreversible_reactions_id, i)
       else
           append!(reversible_reactions_id, i)
       end
   end
   return irreversible_reactions_id, reversible_reactions_id
end

irreversible_reactions_id, reversible_reactions_id = reversibility(myModel)
n_irr = length(irreversible_reactions_id)
n_rev = length(reversible_reactions_id)

function setM(x)
    global M = x
    return
end

function homogenization(lb,ub)
    n = length(lb)
    # Set a large number for M
    setM(+1000.0)
    for i in 1:n
        if lb[i] > 0
            lb[i] = 0
        end
        if ub[i] > 0
            ub[i] = M
        end
        if lb[i] < 0
            lb[i] = -M
        end
        if ub[i] < 0
            ub[i] = 0
        end
    end
    return lb,ub
end

lb, ub = homogenization(lb, ub)

# The Stoichiometric Matrix :

S = Matrix(S)
row_num, col_num = size(S)
println("The Stoichiometric Matrix : $row_num * $col_num")
println("Number of Irreversible Reactions: $n_irr")
println("Number of Reversible Reactions: $n_rev")

# Irreversible Blocked :

lb_u = zeros(n_irr)
ub_u = ones(n_irr)
model = Model(GLPK.Optimizer)
@variable(model, lb[i] <= V[i = 1:n] <= ub[i])
@variable(model, lb_u[i] <= u[i = 1:n_irr] <= ub_u[i])
@constraint(model, S * V .== 0)
objective_function = ones(n_irr)'*u
@objective(model, Max, objective_function)
@constraint(model, [i in 1:n_irr], u[i] <= V[irreversible_reactions_id[i]])
optimize!(model)

blocked_reactions = []
blocked_reactions_index = []
for i in range(1,n_irr)
    if isapprox(value(u[i]), 0.0, atol=1e-6)
        append!(blocked_reactions, value(u[i]))
        append!(blocked_reactions_index, irreversible_reactions_id[i])
    end
end

irr_blocked_num = length(blocked_reactions)

println("Number of Irreversible Blocked reactions: $irr_blocked_num")
println("Irreversible Blocked Index: $blocked_reactions_index")

# S_Transpose :

S_transpose = S'
S_transpose = Matrix(S_transpose)
row_num_trans, col_num_trans = size(S_transpose)
println("transpose : $row_num_trans * $col_num_trans")

# add unit_vectors for reversible reactions :

unit_vector(i,n) = [zeros(i-1); 1 ; zeros(n-i)]

for i in reversible_reactions_id
        a = unit_vector(i, row_num_trans)
        S_transpose = hcat(S_transpose, a)
end

S_transpose_unitVectors = copy(S_transpose)
row_num_trans_unitVectors, col_num_trans_unitVectors = size(S_transpose_unitVectors)
println("The Stoichiometric Matrix + unitVectors : $row_num_trans_unitVectors * $col_num_trans_unitVectors")

# Remove Irrevesible blocked from Stoichiometric Matrix :

S_transpose_unitVectors_noIrrBlocked = []

S_transpose_unitVectors_noIrrBlocked = S_transpose_unitVectors[setdiff(1:end, blocked_reactions_index), :]

S_transpose_unitVectors_noIrrBlocked = Matrix(S_transpose_unitVectors_noIrrBlocked)
row_num_trans_unitVectors_noIrrBlocked, col_num_trans_unitVectors_noIrrBlocked =
                                            size(S_transpose_unitVectors_noIrrBlocked)

println("The Stoichiometric Matrix + unitVectors - Irrevesible blocked  : $row_num_trans_unitVectors_noIrrBlocked * $col_num_trans_unitVectors_noIrrBlocked")

Q, R = qr(S_transpose_unitVectors_noIrrBlocked)
Q = Matrix(Q)
row_r, col_r = size(R)
row_q, col_q = size(Q)
println("S : $row_num_trans_unitVectors_noIrrBlocked * $col_num_trans_unitVectors_noIrrBlocked")
println("Q : $row_q * $col_q")
println("R : $row_r * $col_r")

St_endCol = m
epsilon_row = 10 ^ -15
St = R[:,range(1,St_endCol)]
row_St, col_St = size(St)
println("St : $row_St * $col_St")
global r = 0
for row in eachrow(St)
    r = r + 1
    if norm(row) <= epsilon_row
        break
    end
end
println(r)

R_down = R[range(r,end),range(St_endCol+1,end)]
row_Rdown, col_Rdown = size(R_down)
println("R_down : $row_Rdown * $col_Rdown")

global c = 0
epsilon_col = 10 ^ -15
rev_blocked_col = []
for col in eachcol(R_down)
    c = c + 1
    if norm(col) <= epsilon_col
        append!(rev_blocked_col, c)
    end
end

rev_blocked = []
for i in rev_blocked_col
    append!(rev_blocked, reversible_reactions_id[i])
end

length(rev_blocked) + irr_blocked_num

