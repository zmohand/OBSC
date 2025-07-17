using BenchmarkTools
using Random


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
function is_curtailment_possible(; t_max::Int64, first_curtailment::Int64, last_curtailment::Int64, delta_min::Int64, delta_max::Int64, t_c_b::Int64)
    # Je separe chaque cas où la reduction n'est pas possible (c'est plus lisible)

    if t_c_b == -1 
        return false
    end

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
function init_curtailments(; t_max::Int64, delta::Int64, delta_min::Int64, delta_max::Int64, discharge_levels::Vector{Float64}, p_b::Float64, p_max::Float64, w::Array{Float64})
    vertex_array = Array{Vertex}(undef, 0)
    sequence_id = 1

    # Parcours tous les curtailments possibles (avec niveau de décharge) et les rajoute dans l'ensemble de sommets
    for first_c = 1:t_max
        for last_c = first_c:t_max
            for d_lvl in discharge_levels
                
                t_c_b = calcul_t_c_b(t_max = t_max, c= Vertex(sequence_id, first_c, last_c, d_lvl), delta=delta, p_b = p_b, p_max = p_max, w = w)
                
                if is_curtailment_possible(t_max = t_max, first_curtailment = first_c, last_curtailment = last_c, delta_min = delta_min, delta_max = delta_max, t_c_b = t_c_b)
                    push!(vertex_array, Vertex(sequence_id, first_c, last_c, d_lvl))
                    sequence_id = sequence_id + 1
                end
            end
        
        end
    end

    return vertex_array

end


# Fonction qui va calculer les niveau de d possibles
function first_part_values_d(dmin::Float64, alpha::Float64, b_max::Float64, b_min::Float64, d_min::Float64)
    discharges = Array{Float64}(undef, 0)
    p_set = [i for i in 0:(floor(Int, 1/(1/alpha)))]
    for p in p_set
        push!(discharges, (dmin + p*(1/alpha)*(b_max-b_min-d_min)))
    end
    
    return discharges

end


# Fonction qui va initialiser les sommets en lien avec l'heuristique améliorée
function init_curtailments_improved_heuristic(; t_max::Int64, delta::Int64, delta_min::Int64, delta_max::Int64, alpha::Float64, p_b::Float64, p_max::Float64, w::Array{Float64}, discharge_min::Float64, discharge_max::Float64, p_TSO::Float64, b_max::Float64, b_min::Float64)
    sequence_id = 1
    temp_arrays = Tuple{Int, Int, Float64}[]

    # Tout ça va s'occuper de charger les vertex liés à la premiere partie du set D
    for first_c = 1:t_max
        for last_c = first_c:t_max

            p_c_max = calcul_p_c_j_max( c1 = Vertex(-20, 0, 0, 0), c2 = Vertex(-20, first_c, last_c, 0), delta= delta, p_b = p_b, p_max = p_max, w=w, p_TSO=p_TSO, t_max=t_max, dummy=true)

            dmin = d_min(curtailment = Vertex(-20, first_c, last_c, 0), delta=delta, w=w, p_c_max=p_c_max, discharge_min = discharge_min)

            for d in first_part_values_d(dmin, alpha, b_max, b_min, discharge_min)
                t_c_b = calcul_t_c_b(t_max = t_max, c= Vertex(sequence_id, first_c, last_c, d), delta=delta, p_b = p_b, p_max = p_max, w = w)
                
                if is_curtailment_possible(t_max = t_max, first_curtailment = first_c, last_curtailment = last_c, delta_min = delta_min, delta_max = delta_max, t_c_b = t_c_b)
                    push!(temp_arrays, (first_c, last_c, d))
                    
                end

            end
        end
    end
    # Ajout de la deuxieme partie
    for v in temp_arrays
        dmax = d_max(curtailment=Vertex(-20, v[1], v[2], v[3]), delta=delta, w=w, discharge_max=discharge_max, b_max=b_max, b_min=b_min)
        if !((v[1], v[2], dmax) in temp_arrays)
            t_c_b = calcul_t_c_b(t_max = t_max, c= Vertex(sequence_id, v[1], v[2], dmax), delta=delta, p_b = p_b, p_max = p_max, w = w)
            if is_curtailment_possible(t_max = t_max, first_curtailment = v[1], last_curtailment = v[2], delta_min = delta_min, delta_max = delta_max, t_c_b = t_c_b)
                push!(temp_arrays, (v[1], v[2], dmax)) 
            end
        end
    end

    # Ajout de la 3eme partie
    for v in temp_arrays
        tau_B = v[1] + 1
        while (tau_B <= t_max)
            d = sum(min(p_b, p_max - w[i]) for i in (v[1]+1):tau_B)
            if ! ((v[1], v[2], d) in temp_arrays)
                push!(temp_arrays, (v[1], v[2], d))
            end
            tau_B = tau_B + 1;
        end

    end

    vertex_array = Array{Vertex}(undef, 0)
    # Et on met tout sous la forme d'un array
    for elm in temp_arrays
        push!(vertex_array, Vertex(sequence_id, elm[1], elm[2], elm[3]))
        sequence_id = sequence_id + 1
    end
    return vertex_array 

end


#fonction qui calcule tous les levels discharge possibles
function load_discharge_levels(; b_max::Float64, b_min::Float64, discharge_precision::Int64)

    discharge_array = Array{Float64}(undef, 0)
    level = discharge_precision/100

    while (level < 1)
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
    if t_c_b > t_max
        return -1
    end
    reference_value = delta * sum(min(p_b, p_max - w[t]) for t in (c.curtailment_end+1):t_c_b)

    while reference_value < c.discharge
        t_c_b = t_c_b + 1
        # On vérifie qu'on ne sort pas du time horizon
        if t_c_b > t_max
            return -1
        end
        reference_value = delta * sum(min(p_b, p_max - w[t]) for t in (c.curtailment_end+1):t_c_b)
    end

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
    if (c2.curtailment_start == 1)
        u_f_c_minus_1 = 0
    else
        if dummy
            u_f_c_minus_1 = w[c2.curtailment_start-1]
        else
            u_f_c_minus_1 = w[c2.curtailment_start-1] + calcul_u_b_f_c_minus_one(c1 = c1, c2 = c2, delta = delta, p_b = p_b, p_max = p_max, w = w, t_max = t_max)
        end
    end
    omega_c = (sum(w[t] for t in c2.curtailment_start:c2.curtailment_end) + u_f_c_minus_1) / (c2.curtailment_end - c2.curtailment_start + 2)

    return max(0, omega_c - p_TSO)

end


# Fonction qui calcule la borne sup et la borne inf de d_c_j (qui commence après d_c_i)
function calcul_bornes_d_c_j(; c1::Vertex, c2::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, b_max::Float64, b_min::Float64, t_max::Int64, dummy::Bool)
    p_c_j_max = calcul_p_c_j_max(c1 = c1, c2 = c2, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, t_max = t_max, dummy = dummy)
    
    borne_inf = sum(delta * max(min(w[t], discharge_min), w[t] - p_c_j_max) for t in c2.curtailment_start:c2.curtailment_end)
    
    
    borne_sup = min(
        sum(delta * min(w[t], discharge_max) for t in c2.curtailment_start:c2.curtailment_end),
        b_max - b_min
    )
    
    

    return borne_inf, borne_sup
end


# Fonction qui va vérifier si un arc entre deux sommets est possible
function is_arc_possible(; t_max::Int64, c1::Vertex, c2::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, b_max::Float64, b_min::Float64, dummy::Bool)
    
    if (!dummy)
        t_c1_b = calcul_t_c_b(t_max = t_max, c = c1, delta = delta, p_b = p_b, p_max = p_max, w = w)
        # Cas où on sort des bornes
        if t_c1_b == -1
            return false
        end

        # Cas où la deuxieme reduction commence avant que la batterie ne se recharge
        if (t_c1_b >= c2.curtailment_start)
            return false
        end
    end


    # Il faut maintenant calculer les bornes sup et inf pour d_c_2
   
    borne_inf, borne_sup = calcul_bornes_d_c_j(c1 = c1, c2 = c2, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, b_max = b_max, b_min = b_min, t_max = t_max, dummy=dummy)
    
    if c2.discharge < borne_inf || c2.discharge > borne_sup
        return false
    end

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

    
    return energy_cost[t_c_b]*curtailment.discharge/delta + somme
end


# Fonction qui calcule G'(t) MODE OTR(TODO), ref page 982 pour changer en mode FTR
function calcul_g_prime(; energy_cost::Array{Float64}, reward::Array{Float64}, delta::Int64, time::Int64, OTR::Bool=true)
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

# Fonction qui va calculer d_max
function d_max(; curtailment::Vertex, delta::Int64, w::Array{Float64}, discharge_max::Float64, b_max::Float64, b_min::Float64)
    return min(sum(d_t_max(delta=delta, time=i, w=w, discharge_max=discharge_max) for i in curtailment.curtailment_start:curtailment.curtailment_end) , b_max-b_min)
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
    
  
    g_prime = calcul_g_prime(energy_cost = energy_cost, reward = reward, delta = delta, time = l_set[i_star])
    d_min_i_star = d_min(curtailment = curtailment, delta = delta, w = w, discharge_min = discharge_min, p_c_max = p_c_max)
    somme_du_terme_i_star = sum(begin
        d_t_max(delta = delta, time = l_set[i], w = w, discharge_max = discharge_max) - d_t_min(delta = delta, time = l_set[i], w = w, discharge_min = discharge_min, p_c_max = p_c_max)
    end for i = 1:(i_star-1); init=0)
   
    val_i_star = val_i_star - (g_prime*(d_min_i_star + somme_du_terme_i_star))

    curr_sum = curr_sum + val_i_star

    #Calcul des indices supérieurs (strictement) à i_star
   
    for i = i_star:length(l_set)
        curr_sum = curr_sum + calcul_g_prime(energy_cost = energy_cost, reward = reward, delta = delta, time = l_set[i])*d_t_min(delta = delta, time = l_set[i], w = w, discharge_min = discharge_min, p_c_max = p_c_max)

    end

    return curr_sum
end


# Fonction qui va calculer le poids des arcs (Fonction cout)
function calcul_weight_arcs(; energy_cost::Array{Float64}, reward::Array{Float64}, w::Array{Float64}, t_max::Int64, delta::Int64, curtailment::Vertex,preced_curtailment::Vertex, p_b::Float64, p_max::Float64, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, dummy::Bool)
    battery_charge_cost = calcul_battery_charge_cost(energy_cost = energy_cost, t_max = t_max, delta = delta, curtailment = curtailment, p_b = p_b, p_max = p_max, w = w)
    saving_from_curtailment = calcul_saving_from_curtailment(curtailment = curtailment, preced_curtailment = preced_curtailment, energy_cost = energy_cost, reward = reward, w = w, delta = delta, p_b = p_b, p_max = p_max, p_TSO = p_TSO, t_max = t_max, discharge_min = discharge_min, discharge_max = discharge_max, dummy = dummy)
    return saving_from_curtailment - battery_charge_cost
end


# Fonction qui va renvoyer l'ensemble des arcs du graphe
function init_arcs(; energy_cost::Array{Float64} ,vertex_array::Array{Vertex}, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, b_max::Float64, b_min::Float64, t_max::Int64, reward::Array{Float64})
    n = length(vertex_array)
    arcs_array = [Vector{Arc}() for _ in 1:n]
    # On regarde tous les sommets, si on peut mettre un arc entre 2 sommets on le met
    for s in vertex_array
        for s_prime in vertex_array
            if is_arc_possible(t_max = t_max, c1 = s, c2 = s_prime, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, b_max = b_max, b_min = b_min, dummy = false)
   
                weight = calcul_weight_arcs(energy_cost = energy_cost, reward = reward, w = w, t_max = t_max, delta = delta, curtailment = s_prime, preced_curtailment = s, p_b = p_b, p_max = p_max, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, dummy = false)
                
                push!(arcs_array[s.id], Arc(s, s_prime, weight)) 
            end 
        end
    end

    return arcs_array
end

# Function to add the dummy curtailments
function dummy_curtailment(; vertex_array::Array{Vertex}, arcs_array::Array{Array{Arc, 1}}, t_max::Int64, energy_cost::Array{Float64} , delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, reward::Array{Float64}, b_max::Float64, b_min::Float64)
    push!(vertex_array, Vertex(length(vertex_array)+1, -1, -1, 0))
    push!(vertex_array, Vertex(length(vertex_array) + 1, t_max + 1, t_max + 1, 0))
    append!(arcs_array, [Vector{Arc}() for _ in 1:2])

    for i = 1:(length(vertex_array)-2)
        # Pour v_s (premier sommet du graphe)
        if is_arc_possible(t_max = t_max, c1 = vertex_array[length(vertex_array)-1], c2 = vertex_array[i], delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, b_max = b_max, b_min = b_min, dummy = true)
            weight = calcul_weight_arcs(energy_cost = energy_cost, reward = reward, w = w, t_max = t_max, delta = delta, curtailment = vertex_array[i], preced_curtailment = vertex_array[length(vertex_array)-1], p_b = p_b, p_max = p_max, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, dummy = true)
            push!(arcs_array[length(vertex_array)-1], Arc(vertex_array[length(vertex_array)-1], vertex_array[i], weight))
        end

        # Pour v_t (dernier sommet du graphe)
        push!(arcs_array[length(vertex_array)], Arc(vertex_array[i], vertex_array[length(vertex_array)], 0))
    end

    # On rajoute de v_s à v_t
    push!(arcs_array[length(vertex_array) - 1], Arc(vertex_array[length(vertex_array) - 1], vertex_array[length(vertex_array)], 0))

    start_vertex_id = length(vertex_array) - 1
    end_vertex_id = length(vertex_array)

    return start_vertex_id, end_vertex_id
end

# Fonction qui va calculer la valeur de référence
function calcul_ref_value(w::Array{Float64}, energy_cost::Array{Float64})
    ref_value = 0
    for i = 1:(length(w))
        ref_value = ref_value + w[i]*energy_cost[i]
    end

    return ref_value
end


function init_graph(; improved_heuristic::Bool, t_max::Int64, delta::Int64, delta_min::Int64, delta_max::Int64, discharge_precision::Int64, discharge_min::Float64, discharge_max::Float64, p_TSO::Float64, p_b::Float64, p_max::Float64, b_max::Float64, b_min::Float64, w::Array{Float64}, energy_cost::Array{Float64}, reward::Array{Float64})

    if !improved_heuristic
        discharge_levels_array = load_discharge_levels(b_max = b_max, b_min = b_min, discharge_precision = discharge_precision)
        vertex_array = init_curtailments(t_max = t_max, delta = delta, delta_min = delta_min, delta_max = delta_max, discharge_levels = discharge_levels_array, p_b = p_b, p_max = p_max, w = w)
    else
        vertex_array = init_curtailments_improved_heuristic(t_max=t_max, delta=delta, delta_min = delta_min, delta_max = delta_max, alpha = Float64(discharge_precision), p_b = p_b, p_max = p_max, w= w, discharge_min = discharge_min, discharge_max = discharge_max, p_TSO=p_TSO, b_max = b_max, b_min = b_min)
    end
    arcs_array = init_arcs(energy_cost = energy_cost, vertex_array = vertex_array, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, b_max = b_max, b_min = b_min, t_max = t_max, reward = reward)
    number_array_without_dummy = sum(length(arcs_array[i]) for i in 1:length(vertex_array))
    # On ajoute les dummy curtailment et leurs arcs
    start_v_id, end_v_id = dummy_curtailment(vertex_array = vertex_array, arcs_array = arcs_array, t_max = t_max, energy_cost = energy_cost, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, reward = reward, b_max = b_max, b_min = b_min)
    println("Nombre de sommets avec improved heuristic ? ", improved_heuristic, " nb égal =" ,length(vertex_array))
    # valeur de référence (le cout de l'energie sans réductions)
    reference_value = calcul_ref_value(w, energy_cost)

    return vertex_array, arcs_array, start_v_id, end_v_id, reference_value, number_array_without_dummy
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

function solve(; t_max::Int64, delta::Int64, delta_min::Int64, delta_max::Int64, discharge_precision::Int64, discharge_min::Float64, discharge_max::Float64, p_TSO::Float64, p_b::Float64, p_max::Float64, b_max::Float64, b_min::Float64, w::Array{Float64}, reward::Array{Float64}, energy_cost::Array{Float64})
    v_array, a_array, start_id, end_id, ref_value = init_graph(; t_max=t_max, delta=delta, delta_min = delta_min, delta_max = delta_max, discharge_precision = discharge_precision, discharge_min = discharge_min, discharge_max = discharge_max, p_TSO = p_TSO, p_b = p_b, p_max = p_max, b_max= b_max, b_min = b_min, w = w
, energy_cost = energy_cost, 
reward = reward)
    dist, pred = longest_path_dag(v_array, a_array, start_id)

    return dist, pred, start_id, end_id, v_array, ref_value
end


function run_heuristic(; improved_heuristic ::Bool, t_max::Int64, delta::Int64, delta_min::Int64, delta_max::Int64, discharge_precision::Int64, discharge_min::Float64, discharge_max::Float64, p_TSO::Float64, p_b::Float64, p_max::Float64, b_max::Float64, b_min::Float64, w::Array{Float64}, reward::Array{Float64}, energy_cost::Array{Float64})
    run_time = @elapsed begin
        v_array, a_array, start_id, end_id, ref_value, _ = init_graph(
            improved_heuristic = improved_heuristic, t_max = t_max, delta = delta, delta_min = delta_min, delta_max = delta_max,
            discharge_precision = discharge_precision, discharge_min = discharge_min,
            discharge_max = discharge_max, p_TSO = p_TSO, p_b = p_b, p_max = p_max,
            b_max = b_max, b_min = b_min, w = w, energy_cost = energy_cost, reward = reward
        )
        dist, pred = longest_path_dag(v_array, a_array, start_id)
    end

    
    println("Plus long chemin dans ce graphe : ", dist[end_id])
    println("Temps d'execution : ", run_time)
    println("Ancien cout : ", ref_value, " Nouveau cout : ", (ref_value - dist[end_id]))
    println("Chemin pour l'avoir :")
    nodes_list = trace_path(pred, start_id, end_id, v_array)
    for elm in nodes_list
        println("Id : ", elm.id, " Début reduction : ", elm.curtailment_start, " Fin reduction :", elm.curtailment_end, " Déchargement associé : ", elm.discharge)
    end
end


run_heuristic(improved_heuristic = false, t_max=96, delta=15, delta_min = 1, delta_max = 2, discharge_precision=1, discharge_min = 1.45, discharge_max = 14.5, p_TSO = 5.315, p_b = 6.38, p_max = 31.89, b_max= 106.33, b_min = 53.165, w = [10.507996761717218, 9.733542488512132, 9.662331097944357, 7.97224326726, 6.606087922206213, 8.773672329548887, 9.60626010127904, 9.856003251793606, 9.248740881102975, 11.891793110389175, 10.391533646213217, 11.719437014221056, 8.704960677803792, 8.259145174347815, 9.315124481540888, 13.386755716468052, 9.060394752673709, 10.274518685874014, 12.211452668282316, 6.283198789335201, 10.820185858535353, 9.035532967414596, 8.326332297074845, 7.87520450129465, 11.195546248582295, 8.218228652354782, 12.214054853535025, 8.726102904547798, 7.108264034246719, 9.793387963314384, 11.70109361597738, 10.578641167987449, 11.39126914248145, 9.22644853762128, 11.452580244952715, 12.48640157047186, 10.325215575509661, 9.951981685705116, 10.56198771276182, 6.564360098157062, 10.102874543969666, 10.44963168956906, 9.408717758165354, 8.426025469642054, 11.942124136270886, 13.312167082240453, 9.987508466760662, 11.09800485197126, 9.52859952921557, 10.386142062623094, 11.933153408591394, 12.059946132065233, 11.884699326055664, 12.655700372150378, 10.235484544114358, 9.753526967761127, 12.78708316652121, 12.38331089511357, 10.679915341365021, 11.31319982590546, 8.213565124889476, 13.233791186047466, 9.170882914724071, 10.142717819441764, 11.543630717695233, 12.245091474286136, 8.893397046327896, 10.701467922002182, 13.203418025737102, 13.973763193425679, 10.022992745856564, 13.756808473836678, 13.086674139681627, 13.536322882571167, 12.091015027481081, 11.492379291690263, 11.183774681733151, 10.573965562980868, 11.563461353316134, 10.495738680262052, 8.578821278547759, 16.381005368301874, 10.696985265862413, 10.497089044277876, 10.85843372290127, 13.261612987014699, 7.281136579962315, 8.842239548371019, 7.69322599361125, 13.666772837748502, 12.987821449192515, 14.121461804753514, 10.08520957480766, 15.137866569199335, 11.21839279189488, 8.701238898455218]
, energy_cost = [66.96830659209117, 131.10637569940397, 95.26788754405236, 117.29566687642398, 61.00506655920289, 116.17349709569898, 103.4267715051522, 65.39625326925619, 45.92657274123862, 84.15267902113027, 52.377535192378424, 31.2425382930453, 6.085287370699422, 69.18376673888685, 68.20553235062958, 152.8219881981726, 2.596164479807565, 60.93610530998809, 163.7629337834969, 86.33047055069753, 191.95456602280734, 90.93880490817408, 180.19880964385854, 177.6844727394514, 116.48648767998692, 28.608933560344, 5.612476316465644, 85.88119050131469, 127.09821023849398, 49.49732827414006, 54.73774528966091, 172.8008883419809, 50.091443980246495, 9.824808618836615, 172.98846493089158, 174.1546462864453, 139.64758205664276, 123.4648819667046, 143.3678943325128, 147.75229531411216, 115.05864095315214, 141.8643037353263, 119.47074922691507, 49.70500299079833, 93.18359431327087, 18.860177607721294, 174.04459594075377, 44.406812801139004, 35.83172454134235, 197.29585633737264, 24.727130729213215, 121.59406710597567, 187.13642429997736, 198.00034283627505, 80.56355529531028, 179.12235738529824, 140.30899003218983, 185.17217317876487, 40.17647896072025, 65.89270788930868, 190.4224601590843, 125.84249895965436, 82.88114112251127, 35.58698975522131, 179.48189651418616, 181.1854219627286, 1.502590024661525, 75.20679604274001, 120.55085270824515, 57.077471319790675, 196.26241979925243, 188.29275781832266, 199.00645284192706, 107.85579724514118, 191.99083012264978, 180.99784066569126, 82.93106056593794, 122.00813267411068, 169.97830824408004, 86.63419241704787, 21.749958704401323, 93.94623407472564, 81.24860456578307, 68.40211451706847, 27.32122170425805, 125.76708925601481, 189.52245102480663, 181.26130722459035, 41.62151096460662, 16.496815246199645, 194.86483976694035, 123.8789616061243, 27.18028427662953, 104.41419288027453, 166.82740116777498, 180.39929550149515], 
reward = [68.19782715099605, 13.213080340189252, 138.94947508699263, 127.72787307068619, 85.6120298065787, 10.106978352630385, 72.27014841057569, 60.4106425203684, 134.1594883862214, 100.29747976487339, 138.0166893134052, 78.95765599840719, 117.31918110848213, 65.36565548006425, 128.4041381531086, 63.441492069042006, 57.735871945376104, 24.545864209843348, 117.40039737632993, 21.170461042205087, 12.307783541284557, 7.50445704526741, 29.616541216301172, 55.09781201864526, 94.66512207653484, 110.575357565956, 48.271263955658476, 9.720675777040794, 103.86476092634099, 73.46706209546234, 76.27141356391415, 58.52095384025741, 10.03599677996438, 2.3494489761896404, 113.63264565941311, 79.59510623134672, 117.39312780238465, 107.72157719199545, 33.24824528819315, 84.14785802045677, 117.45982169590918, 75.56351797366193, 3.8851371521499836, 21.474011466382642, 54.03461798084935, 80.05372924138568, 82.93318165564455, 7.3258818056617825, 38.24834496999896, 4.551487795925668, 45.80068964151658, 121.14774297893655, 102.04721123819151, 19.649730728336515, 127.25233255894342, 36.90861648471872, 51.701933679847706, 36.47892203131234, 16.96192293482337, 34.05355468948644, 78.84277432025074, 71.92453403063305, 73.72679939356648, 16.46012298725128, 67.05975718791353, 96.85586532937867, 76.06247915571817, 84.59025935798881, 36.94948862311884, 66.4285045451389, 19.684415571507355, 43.08399484632591, 121.80480695955029, 101.04641394752306, 127.80610974091825, 96.38042642274947, 31.308216944053324, 122.48143976013517, 102.44098203234688, 99.35495236851585, 30.740728503963396, 106.16217215860134, 29.23168486685862, 51.83181072373858, 57.377614848718274, 71.60233178750724, 89.35736833443818, 43.84451894736371, 61.76177507278357, 53.478362832313515, 134.9982991245862, 52.25686061532778, 105.4195192856677, 26.162900958838808, 28.48942833889613, 108.05259553867432])







