include("graph.jl")

#=
    Dans ce fichier on retrouve la fonction principale qui permet de générer des instances de manière
    pseudo aléatoire
=#



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


#=
Fonction qui permet de générer un tableau de valeur avec une certaine moyenne et un écart type
=#
function tableau_moyenne_controlee(n::Int, m::Float64; ecart_type_ratio::Float64 = 0.1)
    # Génère des valeurs avec une moyenne 0 et un écart-type contrôlé
    sigma = m * ecart_type_ratio
    base = randn(n) * sigma
    
    # Décale pour obtenir une moyenne exacte de m
    base .+= m - mean(base)
    return base
end