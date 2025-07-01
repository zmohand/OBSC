using JuMP
using GLPK


struct curtailment
    start_c :: Int64
    end_c :: Int64
end
# Fonction qui va initiliser l'ensemble des curtailment possibles
function init_curtailments(; time_max::Int64, delta_min::Int64, delta_max::Int64)
    result_list = Vector{curtailment}()

    for f_c in 1:time_max
        for l_c in f_c:time_max
            if ((l_c - f_c + 1) <= delta_max) && ((l_c - f_c + 1) >= delta_min)
                push!(result_list, curtailment(f_c, l_c))
            end

        end
    end

    return result_list

end


# Fonction qui va initiliser l'ensemble des curtailment possibles à l'instant t
function init_c_t(; time_max::Int64, curtailments_set::Vector{curtailment})
    result = [[] for _ in 1:time_max]

    for current_time in 1:time_max
        for elm in curtailments_set
            if (elm.start_c <= current_time) && (elm.end_c >= current_time)
                push!(result[current_time], elm)
            end 

        end
    end

    return result
end


model = Model(GLPK.Optimizer)
t = 5
delta_min = 1
delta_max = 2


w = [11.0, 17.0, 12.0, 11.3, 14.2]
energy_cost = [23.0, 12.0, 13.0, 15.0, 10.0]
reward = [4.0, 5.0, 3.4, 3.9, 5.0]
p_max = 20
p_B = 3.0
p_TSO = 1.0
b_min = 30.0
b_max = 45.0

delta = 60

discharge_max = 5.0
discharge_min = 1.0

possible_curtailments = init_curtailments(time_max = t, delta_min = delta_min, delta_max = delta_max)



possible_curtailments_t = init_c_t(time_max = t, curtailments_set = possible_curtailments)

# les variables et contraintes sur les variables

@variable(model, u_D[1:t] >= 0)
@constraint(model, [i=1:t], u_D[i] <= w[i])

@variable(model, u_B[1:t] >= 0)
@constraint(model, [i=1:t], u_B[i] <= p_B)

@variable(model, b_min <= x[1:(t+1)] <= b_max)

@variable(model, z[1:t], Bin)

@variable(model, lin_side[1:t], Bin)

@variable(model, p_c_max[c ∈ possible_curtailments] >= 0)

@variable(model, y_c[c ∈ possible_curtailments], Bin)

@variable(model, 0 <= lin_ypmax_c[c ∈ possible_curtailments] <= p_max)

@variable(model, lin_sidepmax_c[c ∈ possible_curtailments], Bin)

@variable(model, 0 <= lin_xy[1:t, c ∈ possible_curtailments] <= b_max)



# les contraintes

# (5)
@constraint(model, [i=1:t], z[i] == sum(y_c[j] for j in possible_curtailments_t[i]))

# (6)
@constraint(model, [i=1:t], ( x[i] - x[i+1] ) <= (delta * discharge_max * z[i]))

# (7)
@constraint(model, [i=1:t], ( -x[i] + x[i+1] ) >= ( b_max - b_min )*(1 - z[i]))

# (8)
@constraint(model, [i=1:t], ( x[i] - x[i+1] ) >= ( delta * min(w[i], discharge_min)*z[i] - delta*p_B*(1-z[i])))

# (9) JAI UN DOUTE
@constraint(model, [i=1:t], (b_max * sum(y_c[j] for j in possible_curtailments_t[i] if i==j.start_c )) <= x[i])

# (10)
@constraint(model, [i=1:t], ( x[i+1] - x[i] ) == ( delta * (u_B[i] + u_D[i] - w[i])))

# (11)
@constraint(model, x[1] == x[t+1])
@constraint(model, x[t+1] == b_max)
@constraint(model, x[1] == b_max)

# (34)
@constraint(model, [i=1:t],( delta * u_B[i] ) <= ( b_max - x[i] - z[i]*b_max + sum(lin_xy[i, j] for j in possible_curtailments_t[i]) ))

# (35)
@constraint(model, [i=1:t], u_B[i] <= ((1 - z[i]) * min(p_B, p_max - w[i])))

# (36)
@constraint(model, [i=1:t], (b_max / delta - x[i] / delta - min(p_B, p_max - w[i])) 
    <= ( max(p_max, b_max/delta)*lin_side[i])            
)

# (37)
@constraint(model, [i=1:t], ( min(p_B, p_max - w[i]) - b_max/delta + x[i]/delta)
    <= ( max(p_max, b_max/delta)*(1-lin_side[i]) )
)

# (38)
@constraint(model, [i=1:t], ( delta*u_B[i] )
    >= ( b_max - x[i] - z[i]*b_max + sum(lin_xy[i, j] for j in possible_curtailments_t[i]) - max(delta*p_max, b_max)*lin_side[i] )
)

# (39)
@constraint(model, [i=1:t], u_B[i]
    >= ( (1-z[i])*min(p_B, p_max-w[i]) - max(p_max, b_max/delta)*(1-lin_side[i]) )
)

# (13)
@constraint(model, [i=1:t], ( (w[i] - discharge_max)*z[i] + w[i]*(1-z[i]) )
    <= u_D[i]
)

# (40)
@constraint(model, [i=1:t], ( u_D[i] )
    <= ( w[i]*(1-z[i]) + sum(lin_ypmax_c[j] for j in possible_curtailments_t[i]) )
)

# (41)
@constraint(model,
    [c in possible_curtailments],
    p_c_max[c] >= (
        (sum(w[j] for j in max(1, c.start_c-1):c.end_c)   # évite j = 0
         +  x[c.start_c]/delta
         -  (c.start_c == 1 ? b_max : x[c.start_c-1]/delta))  # évite x[0]
        / (c.end_c - c.start_c + 2) - p_TSO)
)

# (42)
@constraint(model, [c in possible_curtailments], ( ( (sum(w[j] for j in (max(1, c.start_c - 1)):c.end_c) + x[c.start_c]/delta - (c.start_c == 1 ? b_max : x[c.start_c-1]/delta))
        / (c.end_c - c.start_c + 2) ) - p_TSO ) 
        <= ( 2*p_max*lin_sidepmax_c[c] )
)

# (43)
@constraint(model, [c in possible_curtailments], 
    (-((sum(w[j] for j in (max(1, c.start_c - 1)):c.end_c) + x[c.start_c]/delta - (c.start_c == 1 ? b_max : x[c.start_c-1]/delta)) / (c.end_c - c.start_c + 2)) + p_TSO)
    <= (2*p_max*(1-lin_sidepmax_c[c]))
)

# (44)
@constraint(model, [c in possible_curtailments],
    p_c_max[c]
    <= (((sum(w[j] for j in (max(1, c.start_c - 1)):c.end_c) + x[c.start_c]/delta - (c.start_c == 1 ? b_max : x[c.start_c-1]/delta)) / (c.end_c - c.start_c + 2)) -p_TSO + 2*p_max*(1-lin_sidepmax_c[c]))
)

# (45)
@constraint(model, [c in possible_curtailments], p_c_max[c] <= (2*p_max*lin_sidepmax_c[c]))

# (46)
@constraint(model, [c in possible_curtailments, i=1:t], lin_xy[i, c] <= (y_c[c]*b_max))

# (47)
@constraint(model, [c in possible_curtailments, i=1:t], lin_xy[i, c] <= x[i])

# (48)
@constraint(model, [c in possible_curtailments, i=1:t], lin_xy[i, c] 
    >= ( x[i] - (1-y_c[c])*b_max)
)

# (49)
@constraint(model, [c in possible_curtailments], lin_ypmax_c[c] <= (y_c[c]*p_max))

# (50)
@constraint(model, [c in possible_curtailments], lin_ypmax_c[c] <= p_c_max[c])

# (51)
@constraint(model, [c in possible_curtailments], lin_ypmax_c[c] <= (p_c_max[c] - (1 - y_c[c])*p_max))




# 3. Définir la fonction objectif (minimiser ou maximiser)
@objective(model, Min, 
    sum(energy_cost[i]*(u_B[i]+ u_D[i]) for i in 1:t) 
    - sum(reward[j]*(w[j] - u_D[j]) for j in 1:t) #OTR
)

# 5. Résoudre
optimize!(model)

# 6. Vérifier le statut et récupérer les résultats
status = termination_status(model)
if status == MOI.OPTIMAL
    println("Solution optimale trouvée :")
else
    println("Pas de solution optimale (status : $status)")
end
