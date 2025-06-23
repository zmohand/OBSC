struct Vertex
    id :: Int64 
    curtailment_start :: Int64
    curtailment_end :: Int64
    discharge :: Float64
end


struct Arc
    curtailment_A :: Vertex
    curtailment_B :: Vertex
    weight :: Int64
end

struct Graph
    vertex_list :: Array{Vertex}
    arc_list :: Array{Array{Arc, 1}}
end


#= Implémentation de Bellman-Ford

"""function Bellman_Ford(g :: Graph, s :: Vertex)
    distance = fill(Inf64, (length(g.vertex_list), 0))
    
    distance[s.id] = 0

    for k in 1:(length(g.vertex_list) - 1)
        for 
    end


end=#


# Fonction qui va vérifier si une réduction donnée (f_c, l_c) est possible 
function is_curtailment_possible(t_max::Int64, first_curtailment::Int64, last_curtailment::Int64, delta_min::Int64, delta_max::Int64)
    # Je separe chaque cas où la reduction n'est pas possible (c'est plus visible)
    if (first_curtailment > last_curtailment)
        return false
    end

    if (first_curtailment >= t_max || last_curtailment > t_max)
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
function init_curtailments(t_max::Int64, delta::Int64, delta_min::Int64, delta_max::Int64, discharge_levels::Vector{Float64})
    vertex_array = Array{Vertex}(undef, 0)
    sequence_id = 1

    # Parcours tous les curtailments possibles (avec niveau de décharge) et les rajoute dans l'ensemble de sommets
    for first_c = 1:t_max
        for last_c = first_c:t_max
            if is_curtailment_possible(t_max, first_c, last_c, delta_min, delta_max)
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
function load_discharge_levels(b_max::Float64, b_min::Float64, discharge_precision::Int64)

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
function calcul_t_c_b(t_max::Int64, c::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64})
    t_c_b = c.curtailment_end+1
    # On vérifie qu'on ne sort pas du time horizon
    if t_c_b > t_max
        return -1
    end
    reference_value = delta*sum(min(p_b, p_max - w[t]) for t = (c.curtailment_end+1):t_c_b)

    while reference_value < c.discharge
        t_c_b = t_c_b + 1
        # On vérifie qu'on ne sort pas du time horizon
        if t_c_b > t_max
            return -1
        end
        reference_value = delta * sum((min(p_b, p_max - w[t]) for t = (c.curtailment_end+1):t_c_b); init=0.0)
    end

    #println("Valeur retournée par t_c_b = ",t_c_b, " appel sur curtailment commencant à ", c.curtailment_start, " et finissant à ", c.curtailment_end)
    return t_c_b

end

# Fonction qui va calculer r_c (la quantité d'energie achetée à la derniere periode de recharge)
function calcul_r_c(c::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, t_max::Int64)
    t_c_b = calcul_t_c_b(t_max, c, delta, p_b, p_max, w)
    return (c.discharge/delta) - sum(min(p_b, p_max - w[t]) for t = (c.curtailment_end+1):(t_c_b - 1); init= 0.0)

end

#= Fonction qui va calculer le u_b_f_c-1, c'est à dire, combien on a acheté
d'energie pour recharger la batterie avant le début du curtailment c2, 
peut dépendre du curtailment précendent c1
=#
function calcul_u_b_f_c_minus_one(c1::Vertex, c2::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, t_max::Int64)
    t_c1_b = calcul_t_c_b(t_max, c1, delta, p_b, p_max, w)
    
    # Cas où on le c2 ne dépend pas de c1, pas de corrélation
    if (t_c1_b < (c2.curtailment_start - 1))
        return 0
    else
        return calcul_r_c(c1, delta, p_b, p_max, w, t_max)
    end

end



# Fonction qui va calculer le P_c_j_max de la reduction c_j qui vient après c_i en utilisant le lemme 1 et propriete 1
function calcul_p_c_j_max(c1::Vertex, c2::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, t_max::Int64)

    
    u_f_c_minus_1 = w[c2.curtailment_start-1] + calcul_u_b_f_c_minus_one(c1, c2, delta, p_b, p_max, w, t_max)

    omega_c = (sum(w[t] for t = c2.curtailment_start:c2.curtailment_end) + u_f_c_minus_1) / (c2.curtailment_end - c2.curtailment_start + 2)

    #println("Pcj max achetable = ", max(0, omega_c-p_TSO), " avec c1 commence à ", c1.curtailment_start, " et finissant à ", c1.curtailment_end, " et c2 commencant à ", c2.curtailment_start, " et finissant à ", c2.curtailment_end, " et discharge c1: ", c1.discharge, " et celle de c2 :", c2.discharge)
    return max(0, omega_c - p_TSO)

end


# Fonction qui calcule la borne sup et la borne inf de d_c_j (qui commence après d_c_i)
function calcul_bornes_d_c_j(c1::Vertex, c2::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, b_max::Float64, b_min::Float64, t_max::Int64)
    p_c_j_max = calcul_p_c_j_max(c1, c2, delta, p_b, p_max, w, p_TSO, t_max)
    borne_inf = sum((delta*max(min(w[t], discharge_min), w[t] - p_c_j_max)) for t = c2.curtailment_start:c2.curtailment_end)

    borne_sup = min(sum(delta*min(w[t], discharge_max) for t = c2.curtailment_start:c2.curtailment_end), b_max - b_min)

    #println("Bornes inf =", borne_inf, " et borne sup :", borne_sup)
    return borne_inf, borne_sup
end


# Fonction qui va vérifier si un arc entre deux sommets est possible
function is_arc_possible(t_max::Int64, c1::Vertex, c2::Vertex, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, b_max::Float64, b_min::Float64)
   
    t_c1_b = calcul_t_c_b(t_max, c1, delta, p_b, p_max, w)
    # Cas où on sort des bornes
    if t_c1_b == -1
        return false
    end

    # Cas où la deuxieme reduction commence avant que la batterie ne se recharge
    if (t_c1_b >= c2.curtailment_start)
        return false
    end

    # Il faut maintenant calculer les bornes sup et inf pour d_c_2
    borne_inf, borne_sup = calcul_bornes_d_c_j(c1, c2, delta, p_b, p_max, w, p_TSO, discharge_min, discharge_max, b_max, b_min, t_max)
    if c2.discharge < borne_inf || c2.discharge > borne_sup
        return false
    end

    # Si tout est bon, alors l'arc est possible
    return true


end

# Fonction qui calcule le f_c_B (cout de recharge de la batterie)
function calcul_battery_charge_cost(energy_cost::Array{Float64}, t_max::Int64, delta::Int64, curtailment::Vertex, p_b::Float64, p_max::Float64, w::Array{Float64})

    somme = 0
    t_c_b = calcul_t_c_b(t_max, curtailment, delta, p_b, p_max, w)
    for i = (curtailment.curtailment_end+1):(t_c_b - 1)
        somme = somme + (energy_cost[i] - energy_cost[t_c_b])*min(p_b, p_max - w[i])
    end

    return energy_cost[t_c_b]*curtailment.discharge/delta + somme
end


# Fonction qui calcule G'(t) MODE OTR(TODO), ref page 982 pour changer en mode FTR
function calcul_g_prime(energy_cost::Array{Float64}, reward::Array{Float64}, delta::Int64, time::Int64)
    return (energy_cost[time] + reward[time])/delta
end


# Fonction qui va calculer le poids des arcs (Fonction cout)
function calcul_weight_arcs(energy_cost::Array{Float64}, t_max::Int64, delta::Int64, curtailment::Vertex, p_b::Float64, p_max::Float64, w::Array{Float64})
    battery_charge_cost = calcul_battery_charge_cost(energy_cost, t_max, delta, curtailment, p_b, p_max, w)
    
end


# Fonction qui va renvoyer l'ensemble des arcs du graphe
function init_arcs(energy_cost::Array{Float64} ,vertex_array::Array{Vertex}, delta::Int64, p_b::Float64, p_max::Float64, w::Array{Float64}, p_TSO::Float64, discharge_min::Float64, discharge_max::Float64, b_max::Float64, b_min::Float64, t_max::Int64)
    n = length(vertex_array)
    arcs_array = [Vector{Arc}() for _ in 1:n]

    # On regarde tous les sommets, si on peut mettre un arc entre 2 sommets on le met
    for s in vertex_array
        for s_prime in vertex_array
            if is_arc_possible(t_max, s, s_prime, delta, p_b, p_max, w, p_TSO, discharge_min, discharge_max, b_max, b_min)
                # WEIGHT A CHANGER
                push!(arcs_array[s.id], Arc(s, s_prime, calcul_weight_arcs(energy_cost, t_max, delta, s_prime, p_b, p_max, w))) 
            end 
        end
    end

    return arcs_array
end

function init_graph(; t_max::Int64, delta::Int64, delta_min::Int64, delta_max::Int64, discharge_precision::Int64, discharge_min::Float64, discharge_max::Float64, p_TSO::Float64, p_b::Float64, p_max::Float64, w::Array{Float64}, b_max::Float64, b_min::Float64)

    discharge_levels_array = load_discharge_levels(b_max, b_min, discharge_precision)

    vertex_array = init_curtailments(t_max, delta, delta_min, delta_max, discharge_levels_array)
    
    arcs_array = init_arcs(vertex_array, delta, p_b, p_max, w, p_TSO, discharge_min, discharge_max, b_max, b_min, t_max)
    #just for testing
    println(arcs_array)
    return vertex_array, arcs_array
end





open("resultat.txt", "w") do f
    redirect_stdout(f) do
        init_graph(t_max=5, delta=60, delta_min = 1, delta_max = 2, discharge_precision=1, discharge_min = 0.5, discharge_max = 0.8, p_TSO = 0.2, p_b = 2.0, p_max = 3.0, w = [1.0, 0.85, 0.85, 0.85, 1.0], b_max= 180.0, b_min = 60.0)
    end
end