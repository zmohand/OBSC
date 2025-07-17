using JuMP
using GLPK
using MathOptInterface
using HiGHS



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


function solve(t::Int64, delta_min::Int64, delta_max::Int64, delta::Int64, OTR::Bool, b_min::Float64, b_max::Float64, p_TSO::Float64, p_B::Float64, p_max::Float64, energy_cost::Array{Float64}, reward::Array{Float64}, w::Array{Float64})
    model = Model(HiGHS.Optimizer)

    possible_curtailments = init_curtailments(time_max = t, delta_min = delta_min, delta_max = delta_max)

    possible_curtailments_t = init_c_t(time_max = t, curtailments_set = possible_curtailments)

    # les variables et contraintes sur les variables

    @variable(model, u_D[1:t] >= 0)
    @constraint(model, [i=1:t], u_D[i] <= w[i])

    @variable(model, 0 <= u_B[1:t] <= p_B)

    @variable(model, b_min <= x[1:(t+1)] <= b_max)

    @variable(model, z[1:t], Bin)

    @variable(model, lin_side[1:t], Bin)

    @variable(model, p_c_max[c ∈ possible_curtailments] >= 0)

    @variable(model, y_c[c ∈ possible_curtailments], Bin)

    @variable(model, 0 <= lin_ypmax_c[c ∈ possible_curtailments] <= p_max)

    @variable(model, lin_sidepmax_c[c ∈ possible_curtailments], Bin)

    @variable(model, 0 <= lin_xy[1:(t+1), c ∈ possible_curtailments] <= b_max)



    # les contraintes

    # (5)
    @constraint(model, [i=1:t], z[i] == sum(y_c[j] for j in possible_curtailments_t[i]))

    # (6)
    @constraint(model, [i=1:t], ( x[i] - x[i+1] ) <= (delta * discharge_max * z[i]))

    # (7)
    @constraint(model, [i=1:t], ( -x[i] + x[i+1] ) <= ( b_max - b_min )*(1 - z[i]))

    # (8)
    @constraint(model, [i=1:t], ( x[i] - x[i+1] ) >= ( delta * min(w[i], discharge_min)*z[i] - delta*p_B*(1-z[i])))

    # (9) JAI UN DOUTE
    @constraint(model, [i=1:t], (b_max * sum(y_c[j] for j in possible_curtailments_t[i] if i==j.start_c )) <= x[i])

    # (10)
    @constraint(model, [i=1:t], ( x[i+1] - x[i] ) == ( delta * (u_B[i] + u_D[i] - w[i])))

    # (11)
    @constraint(model, x[1] == x[t+1])
    @constraint(model, x[t+1] == b_max)

    # (34)
    @constraint(model, [i=1:t],( delta * u_B[i] ) <= ( b_max - x[i] - z[i]*b_max + sum(lin_xy[i, j] for j in possible_curtailments_t[i]) ))

    # (35)
    @constraint(model, [i=1:t], u_B[i] <= ((1 - z[i]) * min(p_B, p_max - w[i])))

    # (36)
    @constraint(model, [i=1:t], ((b_max / delta) - (x[i] / delta) - min(p_B, p_max - w[i])) 
        <= ( max(p_max, b_max/delta)*lin_side[i])            
    )

    # (37)
    @constraint(model, [i=1:t], ( min(p_B, p_max - w[i]) - (b_max/delta) + (x[i]/delta))
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
            +  (x[c.start_c]/delta)
            -  (c.start_c == 1 ? (b_max/delta) : (x[c.start_c-1]/delta)))  # évite x[0]
            / (c.end_c - c.start_c + 2) - p_TSO)
    )

    # (42)
    @constraint(model, [c in possible_curtailments], ( ( (sum(w[j] for j in (max(1, c.start_c - 1)):c.end_c) + (x[c.start_c]/delta) - (c.start_c == 1 ? (b_max/delta) : (x[c.start_c-1]/delta)))
            / (c.end_c - c.start_c + 2) ) - p_TSO ) 
            <= ( 2*p_max*lin_sidepmax_c[c] )
    )

    # (43)
    @constraint(model, [c in possible_curtailments], 
        (-((sum(w[j] for j in (max(1, c.start_c - 1)):c.end_c) + (x[c.start_c]/delta) - (c.start_c == 1 ? (b_max/delta) : (x[c.start_c-1]/delta))) / (c.end_c - c.start_c + 2)) + p_TSO)
        <= (2*p_max*(1-lin_sidepmax_c[c]))
    )

    # (44)
    @constraint(model, [c in possible_curtailments],
        p_c_max[c]
        <= (((sum(w[j] for j in (max(1, c.start_c - 1)):c.end_c) + (x[c.start_c]/delta) - (c.start_c == 1 ? (b_max/delta) : (x[c.start_c-1]/delta))) / (c.end_c - c.start_c + 2)) -p_TSO + 2*p_max*(1-lin_sidepmax_c[c]))
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
    @constraint(model, [c in possible_curtailments], lin_ypmax_c[c] >= (p_c_max[c] - (1 - y_c[c])*p_max))




    # 3. Fonction objectif

    if (OTR)
        @objective(model, Min, 
            sum(energy_cost[i]*(u_B[i]+ u_D[i]) for i in 1:t) 
            - sum(reward[j]*(w[j] - u_D[j]) for j in 1:t) #OTR
        )
    else
        @objective(model, Min, 
            sum(energy_cost[i]*(u_B[i]+ u_D[i]) for i in 1:t) 
            - sum(((reward[c.start_c]/delta)*(lin_xy[c.start_c, c] - lin_xy[c.end_c+1, c])) for c in possible_curtailments)
        )
    end
    # 5. Résoudre
    run_time = @elapsed begin optimize!(model) end
    
    # 6. Vérifier le statut et récupérer les résultats
    status = termination_status(model)
    primal_status = JuMP.primal_status(model)

    println("Statut de terminaison : ", status)
    println("Statut primal : ", primal_status)

    if primal_status in [JuMP.MOI.FEASIBLE_POINT, JuMP.MOI.NEARLY_FEASIBLE_POINT]
        val = objective_value(model)
        println("Valeur de la fonction objectif : ", val)
        println("Temps d'execution : ", run_time)
        println("Solution (même si sous-optimale) :")
        println("u_D = ", value.(u_D))
        println("u_B = ", value.(u_B))
        println("x   = ", value.(x))
        println("z   = ", value.(z))
        println("y_c = ", Dict(c => value(y_c[c]) for c in possible_curtailments if value(y_c[c]) > 0.5))
        println("p_c_max = ", Dict(c => value(p_c_max[c]) for c in possible_curtailments if value(y_c[c]) > 0.5))
    else
        println("Aucune solution réalisable trouvée.")
    end

end

t = 96
delta = 15
delta_min = 1
delta_max = 2
p_max = 5.46
b_max = 9.11
b_min = 4.555
discharge_min = 0.15
discharge_max = 2.49
p_B = 0.59
p_TSO = 0.91
w = [2.0605216936538686, 2.357993714165484, 0.8786643340750919, 1.7055779798915796, 1.5914480676819962, 2.1186036832390176, 2.497588810069657, 1.9253889611545003, 1.7478711016445072, 1.5630149503371469, 2.510469489053591, 2.4036965496284246, 0.8937965044526882, 2.0770168148605395, 1.938313992236913, 1.3489094746078967, 1.720892551762541, 2.0548365761647256, 2.622681533734383, 2.1222297563754764, 2.1414553541718595, 1.8166852757582401, 1.315811631428889, 1.9354373813725305, 1.7225884387732695, 1.5778167727525982, 1.3707560700877772, 1.861639516112948, 1.908297907971717, 1.5530713436278956, 1.5783902223239452, 1.8147886193042524, 1.6035951175812142, 1.9805915315603362, 1.7268969643963832, 1.3098803160231793, 1.4286974102602739, 1.6959505163153057, 1.8558325284477342, 1.60379934498413, 1.9056898401969313, 2.7652251259437532, 2.199039719396268, 1.9028991200182395, 2.0037988347435616, 1.4446554797463618, 1.6670351022364922, 1.9975252757518067, 1.8727582845289485, 1.900200339114957, 2.0918576345396422, 1.465078906610354, 2.249749022532825, 2.3534672050656584, 2.105713347482375, 1.2396971059848085, 1.7943453181732392, 1.7670921693322579, 2.2155413587324104, 0.8797548637713847, 1.791311549224567, 2.3816822631777166, 1.984147137135304, 1.9713645951547023, 1.1562689040115797, 1.2296636694124556, 1.6681021893454233, 1.3143406595431903, 1.9346342534452816, 1.5737098124997875, 1.825788421215909, 2.3906972085786187, 1.681157000143216, 1.8162748605490056, 2.2122749911661157, 2.008156322078934, 1.5514125594508599, 2.5135920642863216, 1.3798384368811678, 2.1860125721135693, 2.1348028964811245, 2.1765132437515, 2.191268666312925, 2.1685829847493547, 2.2846935631816407, 0.997442984245806, 1.5388390338563382, 2.073330299089844, 1.2628622828944738, 1.1544893223770014, 0.999259593786218, 1.7554160724614907, 2.146305160487601, 2.3118192845415155, 1.1935992944400022, 1.9977229919826698]
energy_cost = [121.02941326586429, 151.53212965514604, 80.14541819006995, 96.24895929632403, 79.91263057422567, 136.0340343102437, 171.01041274651112, 120.1307864597865, 77.52624834631357, 36.57309503650068, 171.1646617428792, 190.62085332511973, 43.98996312621394, 100.76194210169302, 165.83898823181516, 63.3200500132517, 166.5789449753441, 90.9205514840113, 140.00241639199763, 110.70002123309749, 107.4560840637453, 78.25934951581586, 69.88252613046119, 20.18453515806009, 99.72704244898951, 73.57515485936952, 182.5926955960181, 172.69246939724042, 49.66015910864178, 110.79236762735259, 28.833884943724314, 26.404730668431956, 101.0074273269061, 42.69826735556226, 166.62531056550105, 53.90716933262874, 190.98368399603416, 37.990846177798765, 11.333318235870639, 89.93599052979926, 80.56424513298066, 36.19576946231487, 103.01231630258658, 39.47845891727644, 141.5632453596932, 194.04451260036555, 128.53763391841653, 113.43403827912205, 20.783984078762863, 41.51120531294435, 165.15952531816504, 50.385361371068775, 40.49656296596956, 27.149023423553057, 46.413058567964235, 94.4370365166404, 145.63369387137877, 182.26982833406902, 18.62302796667944, 169.0298429636608, 45.49834225036418, 157.75004337555777, 132.20689520805314, 52.190518838912745, 162.27436094316525, 43.6415445126145, 159.95038929297257, 143.4924841083733, 95.25707796215019, 44.0567394406623, 47.91464106860398, 125.74323362165484, 73.79766044962095, 24.007731302669146, 126.24770661938192, 84.33972012507172, 185.46210734450185, 28.3881355576776, 95.84520275703389, 70.9242758114604, 190.10859802887455, 50.362333406359014, 104.10247538301071, 190.66887167258415, 53.79231365146035, 142.90034951380753, 53.84060046686905, 127.30791802953658, 187.04696668171877, 172.35824248553155, 104.03108064119935, 153.86291932801424, 181.59578036724184, 57.31119190655699, 15.86635420573545, 18.207048018106683]
reward = [2.867746098183193, 36.47686605848005, 29.331790752090615, 44.68753832260542, 16.06701849734177, 138.04490345160147, 32.27736409676464, 133.42475464915577, 55.671852777640304, 24.15594803054941, 58.819450677136444, 53.47394609837764, 119.69836199480683, 63.51911743573959, 7.4120206182376585, 44.748958974881944, 109.61593816724341, 24.301347485696034, 110.14390468781258, 53.29074687408956, 113.55144879158838, 26.779021176915506, 17.096815501156236, 125.40590961142502, 132.73756731125852, 75.38796011218777, 66.46711267212751, 124.28539172141875, 73.88241973807672, 53.66939744674815, 19.826680946253877, 119.37829571233193, 118.38733501789412, 138.15886493296762, 117.95668944140809, 101.1469454220555, 50.20824141346058, 100.6467322651785, 71.02742794019775, 96.04771689045587, 29.571928313542063, 27.84179893176111, 31.264287719783788, 57.0810065374923, 138.11279947991798, 34.12412544035623, 33.33897928095696, 58.788948079899704, 111.20278553892452, 138.8358663942162, 12.531220422634707, 10.911685641673069, 35.70039972902015, 19.275388597383202, 102.99861245506139, 88.23818810257505, 93.05791825565655, 17.68060476818899, 109.66641236014938, 44.94509714650726, 67.72609463454727, 55.31253093182447, 68.24092254648595, 119.56012232129011, 5.575705869150971, 20.443584960895027, 74.94445641041368, 30.12527431438077, 41.13075104282669, 18.874914899840377, 132.23547106182772, 26.805311194784533, 109.1738696786712, 74.89325758197432, 100.43837711813237, 103.81853455024564, 62.54277040891591, 95.88596338251091, 74.88953479000766, 126.46169966490541, 38.036634515247, 41.20023106427441, 27.967648617052713, 107.85417830601307, 41.85473844358416, 83.72016098600955, 124.48083669590544, 89.1740918801906, 103.11685276652514, 63.26088614360342, 106.9103514782269, 110.90303182579575, 129.46710915584427, 93.90727997306536, 62.55372425682123, 62.27172982168951]

solve(t, delta_min, delta_max, delta, false, b_min, b_max, p_TSO, p_B, p_max, energy_cost, reward, w)

