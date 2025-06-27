struct Vertex
    id :: Int64 
    curtailment_start :: Int64
    curtailment_end :: Int64
    discharge :: Float64
end


struct Arc
    curtailment_A :: Vertex
    curtailment_B :: Vertex
    weight :: Float64
end

struct Graph
    vertex_list :: Array{Vertex}
    arc_list :: Array{Array{Arc, 1}}
end


# Fonction qui va vérifier si une réduction donnée (f_c, l_c) est possible 
function is_curtailment_possible(; t_max::Int64, first_curtailment::Int64, last_curtailment::Int64, delta_min::Int64, delta_max::Int64)
    # Je separe chaque cas où la reduction n'est pas possible (c'est plus visible)
    if (first_curtailment > last_curtailment)
        return false
    end

    if (first_curtailment >= t_max || last_curtailment >= t_max)
        return false
    end

    if (last_curtailment - first_curtailment + 1) > delta_max
        return false
    end

    if (last_curtailment - first_curtailment + 1) < delta_min
        return false
    end

    return true
end


#Fonction qui va génerer l'ensemble des sommets du graphe
function init_curtailments(; t_max::Int64, delta::Int64, delta_min::Int64, delta_max::Int64, discharge_levels::Vector{Float64})
    vertex_array = Array{Vertex}(undef, 0)
    sequence_id = 1

    # Parcours tous les curtailments possibles (avec niveau de décharge) et les rajoute dans l'ensemble de sommets
    for first_c = 1:t_max
        for last_c = first_c:t_max
            if is_curtailment_possible(t_max = t_max, first_curtailment = first_c, last_curtailment = last_c, delta_min = delta_min, delta_max = delta_max)
                for d_lvl in discharge_levels
                    push!(vertex_array, Vertex(sequence_id, first_c, last_c, d_lvl))
                    sequence_id = sequence_id + 1
                end
            end
        end
    end

    return vertex_array

end


#fonction qui calcule tous les levels discharge possibles
function load_discharge_levels(; b_max::Float64, b_min::Float64, discharge_precision::Int64)

    discharge_array = Array{Float64}(undef, 0)
    level = discharge_precision/100

    while (level <= 1)
        curr_lvl = round(level*b_max)

        # Si cette décharge est bien "possible" on la met dans l'array
        if (b_max - curr_lvl) >= b_min
            push!(discharge_array, curr_lvl)
        else
            # Dans ce cas, plus besoin de continuer
            break
        end

        level = level + discharge_precision/100
    end

    return discharge_array

end


# Fonction qui va calculer le t_c_b d'une réduction c en utilisant la "property 1" 
function calcul_t_c_b(; t_max::Int64, c::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64})
    t_c_b = c.curtailment_end+1
    
    reference_value = delta * sum(min(p_b, p_max - w[t]) for t in (c.curtailment_end+1):t_c_b)

    while reference_value < c.discharge
        t_c_b = t_c_b + 1
        # On vérifie qu'on ne sort pas du time horizon
        if t_c_b > t_max
            return -1
        end
        reference_value = delta * sum(min(p_b, p_max - w[t]) for t in (c.curtailment_end+1):t_c_b)
    end

    #println("Valeur retournée par t_c_b = ",t_c_b, " appel sur curtailment commencant à ", c.curtailment_start, " et finissant à ", c.curtailment_end)
    return t_c_b

end

# Fonction qui va calculer r_c (la quantité d'energie achetée à la derniere periode de recharge)
function calcul_r_c(; c::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, t_max::Int64)
    t_c_b = calcul_t_c_b(t_max = t_max, c = c, delta = delta, p_b = p_b, p_max = p_max, w = w)
    return (c.discharge/delta) - sum(min(p_b, p_max - w[t]) for t = (c.curtailment_end+1):(t_c_b - 1); init= 0.0)

end

#= Fonction qui va calculer le u_b_f_c-1, c'est à dire, combien on a acheté
d'energie pour recharger la batterie avant le début du curtailment c2, 
peut dépendre du curtailment précendent c1
=#
function calcul_u_b_f_c_minus_one(; c1::Vertex, c2::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, t_max::Int64)
    t_c1_b = calcul_t_c_b(t_max = t_max, c = c1, delta = delta, p_b = p_b, p_max = p_max, w = w)
    
    # Cas où on le c2 ne dépend pas de c1, pas de corrélation
    if (t_c1_b < (c2.curtailment_start - 1))
        return 0
    else
        return calcul_r_c(c = c1, delta = delta, p_b = p_b, p_max = p_max, w = w, t_max = t_max)
    end

end



# Fonction qui va calculer le P_c_j_max de la reduction c_j qui vient après c_i en utilisant le lemme 1 et propriete 1
function calcul_p_c_j_max(; c1::Vertex, c2::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, t_max::Int64, dummy::Bool)

    # Si on traite un sommet factice :
    if dummy
        u_f_c_minus_1 = 0
    else
        u_f_c_minus_1 = w[c2.curtailment_start-1] + calcul_u_b_f_c_minus_one(c1 = c1, c2 = c2, delta = delta, p_b = p_b, p_max = p_max, w = w, t_max = t_max)
    end
    omega_c = (sum(w[t] for t in c2.curtailment_start:c2.curtailment_end) + u_f_c_minus_1) / (c2.curtailment_end - c2.curtailment_start + 2)

    #println("Pcj max achetable = ", max(0, omega_c-p_TSO), " avec c1 commence à ", c1.curtailment_start, " et finissant à ", c1.curtailment_end, " et c2 commencant à ", c2.curtailment_start, " et finissant à ", c2.curtailment_end, " et discharge c1: ", c1.discharge, " et celle de c2 :", c2.discharge)
    return max(0, omega_c - p_TSO)

end


# Fonction qui calcule la borne sup et la borne inf de d_c_j (qui commence après d_c_i)
function calcul_bornes_d_c_j(; c1::Vertex, c2::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, b_max::Float64, b_min::Float64, t_max::Int64)
    p_c_j_max = calcul_p_c_j_max(c1 = c1, c2 = c2, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, t_max = t_max, dummy = false)
    borne_inf = sum(delta * max(min(w[t], discharge_min), w[t] - p_c_j_max) for t in c2.curtailment_start:c2.curtailment_end)

    borne_sup = min(
        sum(delta * min(w[t], discharge_max) for t in c2.curtailment_start:c2.curtailment_end),
        b_max - b_min
    )

    #println("Bornes inf =", borne_inf, " et borne sup :", borne_sup)
    return borne_inf, borne_sup
end


# Fonction qui va vérifier si un arc entre deux sommets est possible
function is_arc_possible(; t_max::Int64, c1::Vertex, c2::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, b_max::Float64, b_min::Float64)
   
    t_c1_b = calcul_t_c_b( t_max = t_max, c = c1, delta = delta, p_b = p_b, p_max = p_max, w = w)
    # Cas où on sort des bornes
    println("Premier test")
    if t_c1_b == -1
        return false
    end
    println("Deuxieme test")

    # Cas où la deuxieme reduction commence avant que la batterie ne se recharge
    if (t_c1_b >= c2.curtailment_start)
        return false
    end
    println("Troisieme test")

    # Il faut maintenant calculer les bornes sup et inf pour d_c_2
    borne_inf, borne_sup = calcul_bornes_d_c_j(c1 = c1, c2 = c2, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, b_max = b_max, b_min = b_min, t_max = t_max)
    if c2.discharge < borne_inf || c2.discharge > borne_sup
        return false
    end
    println("4eme test")

    # Si tout est bon, alors l'arc est possible
    return true


end

# Fonction qui calcule le f_c_B (cout de recharge de la batterie)
function calcul_battery_charge_cost(; energy_cost::Array{Float64}, t_max::Int64, delta::Int64, curtailment::Vertex, p_b::Float64, p_max::Float64, w::Array{Float64})

    somme = 0
    t_c_b = calcul_t_c_b(t_max = t_max, c = curtailment, delta = delta, p_b = p_b, p_max = p_max, w = w)
    for i in (curtailment.curtailment_end+1):(t_c_b - 1)
        somme = somme + (energy_cost[i] - energy_cost[t_c_b])*min(p_b, p_max - w[i])
    end
    #println("Prob t_c_b ?", t_c_b, " t_max : ", t_max, " c_start :", curtailment.curtailment_start, " c_end : ", curtailment.curtailment_end)
    return energy_cost[t_c_b]*curtailment.discharge/delta + somme
end


# Fonction qui calcule G'(t) MODE OTR(TODO), ref page 982 pour changer en mode FTR
function calcul_g_prime(; energy_cost::Array{Float64}, reward::Array{Float64}, delta::Int64, time::Int64, OTR::Bool=false)
    if OTR
        return (energy_cost[time] + reward[time])/delta
    else
        return 0 # TODO FTR
    end
end

# Fonction qui va calculer l'ensemble L (ref preuve théorème 1 page 982)
function calcul_array_L(; energy_cost::Array{Float64}, reward::Array{Float64}, delta::Int64, curtailment::Vertex)
    L_temp = Vector{Tuple{Float64, Int64}}()
    for i = curtailment.curtailment_start:curtailment.curtailment_end
        push!(L_temp, (calcul_g_prime(energy_cost = energy_cost, reward = reward, delta = delta, time = i), i))
    end
    
    sort!(L_temp; by = first, rev = true)

    return [i for (_, i) in L_temp]

end


# Fonction qui va calculer le d_t_min
function d_t_min(; delta::Int64, time::Int64, w::Array{Float64}, discharge_min::Float64, p_c_max::Float64)
    return delta*max(min(w[time], discharge_min), w[time] - p_c_max)
end

# Fonction qui va calculer le d_t_max
function d_t_max(; delta::Int64, time::Int64, w::Array{Float64}, discharge_max::Float64)
    return delta*min(w[time], discharge_max)
end

# Fonction qui va calculer le d_min (déchargement minimal associé à un curtailment)
function d_min(; curtailment::Vertex, delta::Int64, w::Array{Float64}, discharge_min::Float64, p_c_max::Float64)
    sum = 0
    for i = curtailment.curtailment_start:curtailment.curtailment_end
        sum = sum + d_t_min(delta = delta, time = i, w = w, discharge_min = discharge_min, p_c_max = p_c_max)
    end
    return sum
end

# Fonction qui calculer le i* (ref théorème 1 page 984)
function calcul_i_star(; l_set::Array{Int64}, curtailment::Vertex, preced_curtailment::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, t_max::Int64, discharge_min::Float64, discharge_max::Float64, dummy::Bool)
    index = 1
    p_c_max = calcul_p_c_j_max(c1 = preced_curtailment, c2 = curtailment, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, t_max = t_max, dummy = dummy)
    value = d_min(curtailment = curtailment, delta = delta, w = w, discharge_min = discharge_min, p_c_max = p_c_max) + sum(
        begin
            d_t_max(delta = delta, time = l_set[i], w = w, discharge_max = discharge_max) -
            d_t_min(delta = delta, time = l_set[i], w = w, discharge_min = discharge_min, p_c_max = p_c_max)
        end for i in 1:index)

    while (curtailment.discharge >= value && (index < length(l_set)))
        index = index + 1
        #println("index ? :", index, " c start:", curtailment.curtailment_end, "c end ?", curtailment.curtailment_end)
        value = d_min(curtailment = curtailment, delta = delta, w = w, discharge_min = discharge_min, p_c_max = p_c_max) + sum(
        begin
            d_t_max(delta = delta, time = l_set[i], w = w, discharge_max = discharge_max) -
            d_t_min(delta = delta, time = l_set[i], w = w, discharge_min = discharge_min, p_c_max = p_c_max)
        end for i in 1:index)
    end


    return index
end


function calcul_saving_from_curtailment(; curtailment::Vertex, preced_curtailment::Vertex, energy_cost::Array{Float64}, reward::Array{Float64}, w::Array{Float64}, delta::Int64, p_b::Float64, p_max::Float64, p_TSO::Float64, t_max::Int64, discharge_min::Float64, discharge_max::Float64, dummy::Bool)
    l_set = calcul_array_L(energy_cost = energy_cost, reward = reward, delta = delta, curtailment = curtailment)
    i_star = calcul_i_star(l_set = l_set, curtailment = curtailment, preced_curtailment = preced_curtailment, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, t_max = t_max, discharge_min = discharge_min, discharge_max = discharge_max, dummy = dummy)

    # Calcul des indices inferieurs à i_star
    curr_sum = 0
    for i = 1:i_star-1
        curr_sum = curr_sum + calcul_g_prime(energy_cost = energy_cost, reward = reward, delta = delta, time = l_set[i])*d_t_max(delta = delta, time = l_set[i], w = w, discharge_max = discharge_max)
    end

    p_c_max = calcul_p_c_j_max(c1 = preced_curtailment, c2 = curtailment, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, t_max = t_max, dummy = dummy)
    #Calcul de l'indice i_star
    val_i_star = calcul_g_prime(energy_cost = energy_cost, reward = reward, delta = delta, time = l_set[i_star])*curtailment.discharge
    val_i_star = val_i_star - (calcul_g_prime(energy_cost = energy_cost, reward = reward, delta = delta, time = l_set[i_star]) * ( d_min(curtailment = curtailment, delta = delta, w = w, discharge_min = discharge_min, p_c_max = p_c_max) + sum(begin
        d_t_max(delta = delta, time = l_set[i], w = w, discharge_max = discharge_max) - d_t_min(delta = delta, time = l_set[i], w = w, discharge_min = discharge_min, p_c_max = p_c_max)
    end for i = 1:(i_star-1); init=0)))

    curr_sum = curr_sum + val_i_star

    #Calcul des indices supérieurs (strictement) à i_star

    for i = i_star+1:length(l_set)
        curr_sum = curr_sum + calcul_g_prime(energy_cost = energy_cost, reward = reward, delta = delta, time = l_set[i])*d_t_min(delta = delta, time = l_set[i], w = w, discharge_min = discharge_min, p_c_max = p_c_max)

    end


    return curr_sum
end


# Fonction qui va calculer le poids des arcs (Fonction cout)
function calcul_weight_arcs(; energy_cost::Array{Float64}, reward::Array{Float64}, w::Array{Float64}, t_max::Int64, delta::Int64, curtailment::Vertex,preced_curtailment::Vertex, p_b::Float64, p_max::Float64, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, dummy::Bool)
    battery_charge_cost = calcul_battery_charge_cost(energy_cost = energy_cost, t_max = t_max, delta = delta, curtailment = curtailment, p_b = p_b, p_max = p_max, w = w)
    #println("BATTERY COST CHARGING : ", battery_charge_cost)
    saving_from_curtailment = calcul_saving_from_curtailment(curtailment = curtailment, preced_curtailment = preced_curtailment, energy_cost = energy_cost, reward = reward, w = w, delta = delta, p_b = p_b, p_max = p_max, p_TSO = p_TSO, t_max = t_max, discharge_min = discharge_min, discharge_max = discharge_max, dummy = dummy)
    #println("SAVING FROM C: ", saving_from_curtailment)
    return saving_from_curtailment - battery_charge_cost
end


# Fonction qui va renvoyer l'ensemble des arcs du graphe
function init_arcs(; energy_cost::Array{Float64} ,vertex_array::Array{Vertex}, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, b_max::Float64, b_min::Float64, t_max::Int64, reward::Array{Float64})
    n = length(vertex_array)
    arcs_array = [Vector{Arc}() for _ in 1:n]

    # On regarde tous les sommets, si on peut mettre un arc entre 2 sommets on le met
    for s in vertex_array
        for s_prime in vertex_array
            if is_arc_possible(t_max = t_max, c1 = s, c2 = s_prime, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, b_max = b_max, b_min = b_min)
                # WEIGHT A CHANGER
                weight = calcul_weight_arcs(energy_cost = energy_cost, reward = reward, w = w, t_max = t_max, delta = delta, curtailment = s_prime, preced_curtailment = s, p_b = p_b, p_max = p_max, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, dummy = false)
                println("Weight : ", weight)
                #println("WEIGHT EN SORTIE : ", weight)
                push!(arcs_array[s.id], Arc(s, s_prime, weight)) 
            end 
        end
    end

    return arcs_array
end

# Function to add the dummy curtailments
function dummy_curtailment(; vertex_array::Array{Vertex}, arcs_array::Array{Array{Arc, 1}}, t_max::Int64, energy_cost::Array{Float64} , delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, reward::Array{Float64})
    push!(vertex_array, Vertex(length(vertex_array)+1, -1, -1, 0))
    push!(vertex_array, Vertex(length(vertex_array) + 1, t_max + 1, t_max + 1, 0))
    append!(arcs_array, [Vector{Arc}() for _ in 1:2])

    for i = 1:(length(vertex_array)-2)
        # Pour v_s (premier sommet du graphe)
        weight = calcul_weight_arcs(energy_cost = energy_cost, reward = reward, w = w, t_max = t_max, delta = delta, curtailment = vertex_array[i], preced_curtailment = vertex_array[length(vertex_array)-1], p_b = p_b, p_max = p_max, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, dummy = true)
        push!(arcs_array[length(vertex_array)-1], Arc(vertex_array[length(vertex_array)-1], vertex_array[i], weight))

        # Pour v_t (dernier sommet du graphe)
        push!(arcs_array[length(vertex_array)], Arc(vertex_array[i], vertex_array[length(vertex_array)], 0))
    end

    # On rajoute de v_s à v_t
    push!(arcs_array[length(vertex_array) - 1], Arc(vertex_array[length(vertex_array) - 1], vertex_array[length(vertex_array)], 0))

    start_vertex_id = length(vertex_array) - 1
    end_vertex_id = length(vertex_array)

    return start_vertex_id, end_vertex_id
end


function init_graph(; t_max::Int64, delta::Int64, delta_min::Int64, delta_max::Int64, discharge_precision::Int64, discharge_min::Float64, discharge_max::Float64, p_TSO::Float64, p_b::Float64, p_max::Float64, b_max::Float64, b_min::Float64, w::Array{Float64}, energy_cost::Array{Float64}, reward::Array{Float64})

    discharge_levels_array = load_discharge_levels(b_max = b_max, b_min = b_min, discharge_precision = discharge_precision)

    vertex_array = init_curtailments(t_max = t_max, delta = delta, delta_min = delta_min, delta_max = delta_max, discharge_levels = discharge_levels_array)
    
    arcs_array = init_arcs(energy_cost = energy_cost, vertex_array = vertex_array, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, b_max = b_max, b_min = b_min, t_max = t_max, reward = reward)

    # On ajoute les dummy curtailment et leurs arcs
    start_v_id, end_v_id = dummy_curtailment(vertex_array = vertex_array, arcs_array = arcs_array, t_max = t_max, energy_cost = energy_cost, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, reward = reward)
    #just for testing
    #println(arcs_array)
    return vertex_array, arcs_array, start_v_id, end_v_id
end




# Tri topo + bellman
function longest_path_dag(vertices::Vector{Vertex}, arcs::Vector{Vector{Arc}}, start_id::Int64)

    ids = [v.id for v in vertices]          
    adj = Dict(id => Arc[] for id in ids)         
    indeg = Dict(id => 0    for id in ids)

    for a in arcs
        for elm in a
            push!(adj[elm.curtailment_A.id], elm)
            indeg[elm.curtailment_B.id] += 1
        end
    end

    # tri topologique 
    order = Int[]
    q = [id for id in ids if indeg[id] == 0]
    while !isempty(q)
        u = popfirst!(q)
        push!(order, u)
        for a in adj[u]
            v = a.curtailment_B.id
            indeg[v] -= 1
            indeg[v] == 0 && push!(q, v)
        end
    end

    # Bellman
    dist = Dict(id => -Inf for id in ids)
    pred = Dict(id => 0 for id in ids)
    dist[start_id] = 0.0

    for u in order
        for a in adj[u]
            v = a.curtailment_B.id
            cand = dist[u] + a.weight
            if cand > dist[v]
                dist[v] = cand
                pred[v] = u
            end
        end
    end
    return dist, pred
end

# Fonction qui va retracer le plus long chemin (pour les tests)
function trace_path(pred::Dict{Int64, Int64}, start_id::Int64, final_id::Int64, vertex_list::Array{Vertex})
    result = []
    pushfirst!(result, vertex_list[final_id])
    current = pred[final_id]

    while (current != start_id && current != 0)
        pushfirst!(result, vertex_list[current])
        current = pred[current]
    end


    pushfirst!(result, vertex_list[start_id])


    return result


end


open("resultat.txt", "w") do f
    redirect_stdout(f) do
        v_array, a_array, start_id, end_id = init_graph(t_max=24, delta=60, delta_min = 1, delta_max = 2, discharge_precision=1, discharge_min = 1.45, discharge_max = 14.5, p_TSO = 5.315, p_b = 6.38, p_max = 31.89, b_max= 106.33, b_min = 53.165, w = [12.0, 10.0, 8.4, 7.6, 6.8, 6.4, 6.4, 7.6, 8.4, 9.2, 9.6, 10.0, 10.4, 10.4, 10.8, 11.2, 11.28, 12.0, 12.4, 12.8, 12.8, 14.0, 14.4, 14.0]
, energy_cost = [
    29.32, 26.21, 23.61, 21.04, 21.62, 26.54, 32.89, 38.30,
    34.68, 39.43, 38.31, 37.30, 35.89, 34.41, 32.17, 31.04,
    31.55, 35.29, 38.60, 39.22, 37.34, 34.96, 35.36, 34.65
], 
reward = [
    39.66, 36.30, 30.15, 11.51,  1.07,  0.70, 13.90, 46.51,
    45.30, 59.51, 58.26, 57.00, 51.00, 46.05, 42.61, 40.97,
    46.96, 49.06, 42.67, 41.37, 40.20, 41.15, 43.13, 41.81
])
        dist, pred = longest_path_dag(v_array, a_array, start_id)
        println("Plus long chemin dans ce graphe : ", dist[end_id])
        println("Chemin pour l'avoir :")
        nodes_list = trace_path(pred, start_id, end_id, v_array)
        for elm in nodes_list
            println("Id : ", elm.id, " Début reduction : ", elm.curtailment_start, " Fin reduction :", elm.curtailment_end, " Déchargement associé : ", elm.discharge)
        end

    end
end


