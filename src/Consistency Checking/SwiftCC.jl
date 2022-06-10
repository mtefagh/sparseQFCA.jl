# SwiftCC

using COBREXA
using JuMP
using GLPK
using LinearAlgebra

myModel = load_model(StandardModel,"/Users/iman/Desktop/M.S.Thesis/SourceCode/Data/MetabolicNetworks/e_coli_core.xml")

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

function swiftCC(myModel)
    
    S, Metabolites, Reactions, Genes, m, n, lb, ub = dataOfModel(myModel)
    
    irreversible_reactions_id, reversible_reactions_id = reversibility(myModel)
    n_irr = length(irreversible_reactions_id)
    n_rev = length(reversible_reactions_id)
    
    lb, ub = homogenization(lb, ub) 
    
    # The Stoichiometric Matrix : 

    S = Matrix(S)
    row_num, col_num = size(S)
    
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
    
    irr_blocked_reactions = []
    for i in range(1,n_irr)
        if isapprox(value(u[i]), 0.0, atol=1e-6)
            append!(irr_blocked_reactions, irreversible_reactions_id[i])
        end
    end
    
    irr_blocked_num = length(irr_blocked_reactions)
    
    # Reversible Blocked
    
    # S_Transpose :

    S_transpose = S'
    S_transpose = Matrix(S_transpose)
    row_num_trans, col_num_trans = size(S_transpose)
    
    # Add unit_vectors for reversible reactions : 

    unit_vector(i,n) = [zeros(i-1); 1 ; zeros(n-i)]

    for i in reversible_reactions_id
        a = unit_vector(i, row_num_trans)
        S_transpose = hcat(S_transpose, a)
    end
    S_transpose_unitVectors = copy(S_transpose)
    row_num_trans_unitVectors, col_num_trans_unitVectors = size(S_transpose_unitVectors)
    
    # Remove Irrevesible blocked from Stoichiometric Matrix : 

    S_transpose_unitVectors_noIrrBlocked = []

    S_transpose_unitVectors_noIrrBlocked = S_transpose_unitVectors[setdiff(1:end, irr_blocked_reactions), :]

    S_transpose_unitVectors_noIrrBlocked = Matrix(S_transpose_unitVectors_noIrrBlocked)
    row_num_trans_unitVectors_noIrrBlocked, col_num_trans_unitVectors_noIrrBlocked = 
                                            size(S_transpose_unitVectors_noIrrBlocked)
    
    # Find Reversible blocked
    
    Q, R = qr(S_transpose_unitVectors_noIrrBlocked)
    Q = Matrix(Q)
    row_r, col_r = size(R)
    row_q, col_q = size(Q)
    
    St_endCol = m
    epsilon_row = 10 ^ -6
    St = R[:,range(1,St_endCol)]
    row_St, col_St = size(St)
    
    # Find special row of the S_T 
    
    r = 0
    for row in eachrow(St)
        r = r + 1
        if norm(row)/sqrt(length(row)) <= epsilon_row
            break
        end
    end
    
    R_down = R[range(r,end),range(St_endCol+1,end)]
    row_Rdown, col_Rdown = size(R_down)
    
    # Find specific columns
    
    c = 0

    epsilon_col = 10 ^ -6
    rev_blocked_reactions_col = []
    for col in eachcol(R_down)
        c = c + 1
        if norm(col)/sqrt(length(col)) <= epsilon_col
            append!(rev_blocked_reactions_col, c)
        end
    end

    # Find reversible blocked 
    
    rev_blocked_reactions = []
    for i in rev_blocked_reactions_col
        append!(rev_blocked_reactions, reversible_reactions_id[i])
    end
    
    # Merge reversbile and irreversible reactions
    
    blocked_index = []
    blocked_index = union(rev_blocked_reactions, irr_blocked_reactions)
    
    blocked_names = []
    for i in blocked_index
        r_name = reactions(myModel)[i]
        push!(blocked_names, r_name)
    end
    
    blocked_names = sort(blocked_names)
    return length(irr_blocked_reactions), length(rev_blocked_reactions)
    
end

x, y = swiftCC(myModel)

println("Irreversible Blocked : ")
println(x)
println("Reversible Blocked : ")
println(y)
