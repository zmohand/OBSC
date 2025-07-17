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
function init_curtailments_improved_heuristic(; t_max::Int64, delta::Int64, delta_min::Int64, delta_max::Int64, alpha::float64, p_b::Float64, p_max::Float64, w::Array{Float64}, discharge_min::Float64, discharge_max::Float64)
    sequence_id = 1
    temp_arrays = Tuple{Int, Int, Float64}[]

    # Tout ça va s'occuper de charger les vertex liés à la premiere partie du set D
    for first_c = 1:t_max
        for last_c = first_c:t_max

            p_c_max = calcul_p_c_j_max( c1 = Vertex(-20, 0, 0, 0), c2 = Vertex(-20, first_c, last_c, 0), delta= delta, p_b = p_b, p_max = p_max, w=w, p_TSO=p_TSO, t_max=t_max, dummy=true)

            dmin = d_min(curtailment = Vertex(-20, first_c, last_c, 0), delta=delta, w=w, p_c_max=p_c_max)

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
        d_max = d_max(curtailment=Vertex(-20, v[1], v[2], v[3]), delta=delta, w=w, discharge_max=discharge_max, b_max=b_max, b_min=b_min)
        if !((v[1], v[2], d_max) in temp_arrays)
            t_c_b = calcul_t_c_b(t_max = t_max, c= Vertex(sequence_id, first_c, last_c, d), delta=delta, p_b = p_b, p_max = p_max, w = w)
            if is_curtailment_possible(t_max = t_max, first_curtailment = first_c, last_curtailment = last_c, delta_min = delta_min, delta_max = delta_max, t_c_b = t_c_b)
                push!(temp_arrays, (first_c, last_c, d)) 
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


function init_graph(; t_max::Int64, delta::Int64, delta_min::Int64, delta_max::Int64, discharge_precision::Int64, discharge_min::Float64, discharge_max::Float64, p_TSO::Float64, p_b::Float64, p_max::Float64, b_max::Float64, b_min::Float64, w::Array{Float64}, energy_cost::Array{Float64}, reward::Array{Float64})

    discharge_levels_array = load_discharge_levels(b_max = b_max, b_min = b_min, discharge_precision = discharge_precision)
    vertex_array = init_curtailments(t_max = t_max, delta = delta, delta_min = delta_min, delta_max = delta_max, discharge_levels = discharge_levels_array, p_b = p_b, p_max = p_max, w = w)
    
    arcs_array = init_arcs(energy_cost = energy_cost, vertex_array = vertex_array, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, b_max = b_max, b_min = b_min, t_max = t_max, reward = reward)
    number_array_without_dummy = sum(length(arcs_array[i]) for i in 1:length(vertex_array))
    # On ajoute les dummy curtailment et leurs arcs
    start_v_id, end_v_id = dummy_curtailment(vertex_array = vertex_array, arcs_array = arcs_array, t_max = t_max, energy_cost = energy_cost, delta = delta, p_b = p_b, p_max = p_max, w = w, p_TSO = p_TSO, discharge_min = discharge_min, discharge_max = discharge_max, reward = reward, b_max = b_max, b_min = b_min)

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


function run_heuristic(; t_max::Int64, delta::Int64, delta_min::Int64, delta_max::Int64, discharge_precision::Int64, discharge_min::Float64, discharge_max::Float64, p_TSO::Float64, p_b::Float64, p_max::Float64, b_max::Float64, b_min::Float64, w::Array{Float64}, reward::Array{Float64}, energy_cost::Array{Float64})
    run_time = @elapsed begin
        v_array, a_array, start_id, end_id, ref_value, _ = init_graph(
            t_max = t_max, delta = delta, delta_min = delta_min, delta_max = delta_max,
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


run_heuristic(t_max=5, delta=60, delta_min = 1, delta_max = 2, discharge_precision=1, discharge_min = 1.45, discharge_max = 14.5, p_TSO = 5.315, p_b = 6.38, p_max = 31.89, b_max= 106.33, b_min = 53.165, w = [12.5, 10.0, 8.4, 7.4, 6.0]
, energy_cost = [37.3, 35.2, 26.3, 23.3, 25.3], 
reward = [43.4, 38.2, 30.2, 16.3, 3.5])


# Fonction qui va calculer le nombre de curtailement qui sont liés (corrélation entre 2 curtailments)
function count_correlated_curtailments(; vertices::Vector{Vertex}, arcs::Vector{Vector{Arc}}, t_max::Int64, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, start_id :: Int64, end_id :: Int64)
    counter = 0
    times = []
    for i in 1:length(vertices)
        if i != start_id && i != end_id

            for j in 1:length(arcs[i])
                c1 = arcs[i][j].curtailment_A
                c2 = arcs[i][j].curtailment_B
                t_c_b = calcul_t_c_b(t_max = t_max, c = c1, delta= delta, p_b = p_b , p_max = p_max, w = w)
                if t_c_b == (c2.curtailment_start-1) && !(t_c_b in times)
                    push!(times, t_c_b)
                    counter = counter + 1
                end

            end
        end
    end

    return counter

end

function tableau_moyenne_controlee(n::Int, m::Float64; ecart_type_ratio::Float64 = 0.1)
    # Génère des valeurs avec une moyenne 0 et un écart-type contrôlé
    sigma = m * ecart_type_ratio
    base = randn(n) * sigma
    
    # Décale pour obtenir une moyenne exacte de m
    base .+= m - mean(base)
    return base
end




#= Fonction de génération pseudo-aléatoire
On va considérer les différents sites décrits dans la thèse
Puis générer aléatoirement les w, energy_cost, reward
=#
function random_generation()
    Random.seed!()

    #=
        Description :
        J'ai pris les sites de la thèse, avec un dictionnaire (clé : nom du site, exemple : S1)
        Valeur : un tuple (a, b) avec a une liste des paramètres non array et le b les parametres array
        pour a : [t_max, delta, delta_min, delta_max, p_max, b_max, b_min, discharge_min, discharge_max, p_B, p_TSO]
        pour b : [[w], [energy_cost], [reward]]
    =#

    # Represente le nombre d'instances "satisfaisantes" trouvées

    counter = 0
    w = Float64[]   
    energy_cost = Float64[]
    reward = Float64[]


    

    #TODO rajouter les autres sites
    # Discharge precision mis à 1%
    sites = Dict(
        "S1"=>([96, 15, 1, 2, 11.79, 39.26, 19.63, 0.538, 5.38, 2.34, 1.965], [[], [], []]),
        "S2"=>([96, 15, 1, 2, 1.41, 2.35, 1.175, 0.064, 0.64, 0.14, 0.235], [[], [], []]),
        "S3"=>([96, 15, 1, 2, 5.46, 9.11, 4.555, 0.15, 2.49, 0.59, 0.91], [[], [], []]),
        "S4"=>([96, 15, 1, 2, 31.89, 106.33, 53.165, 1.45, 14.5, 6.38, 5.315], [[], [], []])
    )
    for k in keys(sites)
        counter = 0
        while counter < 5
            w = tableau_moyenne_controlee(Int(sites[k][1][1]), sites[k][1][5]/3; ecart_type_ratio = rand(0.05:0.01:0.2))
            energy_cost = rand(Int(sites[k][1][1])) .* (200.0 - 1.0) .+ 1.0
            reward = rand(Int(sites[k][1][1])) .* (140.0 - 1.0) .+ 1.0
            v_array, a_array, start_id, end_id, ref_value, number_array_without_dummy = init_graph(
                t_max = Int(sites[k][1][1]), delta = Int(sites[k][1][2]), delta_min = Int(sites[k][1][3]), delta_max = Int(sites[k][1][4]),
                discharge_precision = 1, discharge_min = sites[k][1][8],
                discharge_max = sites[k][1][9], p_TSO = sites[k][1][11], p_b = sites[k][1][10], p_max = sites[k][1][5],
                b_max = sites[k][1][6], b_min = sites[k][1][7], w = w, energy_cost = energy_cost, reward = reward
            )
            counter = count_correlated_curtailments(vertices = v_array, arcs=a_array, t_max= Int(sites[k][1][1]), delta = Int(sites[k][1][2]), p_b =sites[k][1][10], p_max = sites[k][1][5], w = w, start_id = start_id, end_id = end_id)
            println("Counter trouvé", counter)
        end
       
        push!(sites[k][2][1], w)
        push!(sites[k][2][2], energy_cost)
        push!(sites[k][2][3], reward)
        println("passed")

    end

    return sites

end


