@time begin
    using MetX
    using JSON
    using CSV
    using DataFrames
    using SparseArrays
    using CairoMakie
    using Clp
    using Ipopt
end

## ..- -- .- . -- -- - . .. .. -. -.- 
#ISSUE
# - the network (at least the pan network) is unstable.
#   - for instance, it is unsensible to ATPM increases.


## ..- -- .- . -- -- - . .. .. -. -.- 
#MARK: Load rxninfo_tsv
let
    paperson = "/Users/Pereiro/.julia/dev/RegulationImageMC_2024/data/PaperSON"
    monkdata = joinpath(paperson, "data", "monkGenomescaleMetabolicNetwork2022")
    # readdir(monkdata)
    rxninfo_json = JSON.parsefile(joinpath(monkdata, "raw.rstb20210236_si_002.json"))
    tsv_file = joinpath(monkdata, rxninfo_json["data"]["Reaction_Info"]["data"]["file"])
    global rxninfo_tsv = CSV.read(tsv_file, DataFrame)
end

## ..- -- .- . -- -- - . .. .. -. -.- 
#MARK: _parse_eq!
function _parse_eq!(rxn_ast::Dict, rxnEq::String)

    rxn_ast["rxnEq"] = rxnEq
    # parse arrow
    lb = nothing
    ub = nothing
    req = ""
    peq = ""
    if contains(rxnEq, "-->")
        req, peq = split(rxnEq, "-->"; keepempty = false)
        lb = 0
        ub = 1000
    elseif contains(rxnEq, "<--")
        req, peq = split(rxnEq, "<--"; keepempty = false)
        lb = -1000
        ub = 0
    elseif contains(rxnEq, "<=>")
        req, peq = split(rxnEq, "<=>"; keepempty = false)
        lb = -1000
        ub = 1000
    end
    rxn_ast["lb"] = lb
    rxn_ast["ub"] = ub
    rxn_ast["stoi"] = Dict{String, Float64}()
    stoi = rxn_ast["stoi"]

    # reactants
    for (dig1, sign) in [
        (strip.(split(req, "+"; keepempty = false)), -1),
        (strip.(split(peq, "+"; keepempty = false)), 1),
    ]
        for pair in dig1
            dig2 = split(pair, " ")
            if length(dig2) == 1
                stoi[dig2[1]] = sign * 1
            end
            if length(dig2) == 2
                stoi[dig2[2]] = sign * parse(Float64, dig2[1])
            end
        end
    end

    return nothing
end

## ..- -- .- . -- -- - . .. .. -. -.- 
#MARK: _parse_rxn_ast
function _parse_rxn_ast(rxninfo_tsv::DataFrame, ridx)
    rxn_ast = Dict{String, Any}()
    
    rxn = rxninfo_tsv[ridx, "Reaction ID"]
    rxn_ast["rxn"] = rxn
    rxnName = rxninfo_tsv[ridx, "Reaction Name"]
    rxn_ast["rxnName"] = rxnName
    subSys = rxninfo_tsv[ridx, "Subsystem"]
    rxn_ast["subSys"] = subSys
    
    rxnEq = rxninfo_tsv[ridx, "Reaction"]
    _parse_eq!(rxn_ast, rxnEq)

    return rxn_ast
end

## ..- -- .- . -- -- - . .. .. -. -.- 
#MARK: pannet_ast
let
    global pannet_ast = Dict{String, Dict}()
    for ridx in 1:size(rxninfo_tsv, 1)
        rxn_ast = _parse_rxn_ast(rxninfo_tsv, ridx)
        pannet_ast[rxn_ast["rxn"]] = rxn_ast
    end
    
    json_file = joinpath(@__DIR__, "pannet_ast.json")
    open(json_file, "w") do io
        JSON.print(io, pannet_ast, 4)
    end
end

## ..- -- .- . -- -- - . .. .. -. -.- 
#MARK: pannet_ast standard COBRA form
let 
    ## ..- -- .- . --
    # indmap
    rxns_idxmap = Dict{String, Int}()
    rxns_idx0 = 0
    mets_idxmap = Dict{String, Int}()
    mets_idx0 = 0
    
    for (rxn, ast) in pannet_ast
        isempty(rxn) && continue
        get!(rxns_idxmap, rxn) do
            rxns_idx0 += 1
            return rxns_idx0
        end
        for (met, s) in ast["stoi"]
            isempty(met) && continue
            get!(mets_idxmap, met) do
                mets_idx0 += 1
                return mets_idx0
            end
        end
    end

    # For adding exchange reactions
    mets_e = filter(keys(mets_idxmap)) do met
        endswith(met, "_e") || return false
        
        rxn = string("EX_", met)
        get!(rxns_idxmap, rxn) do
            rxns_idx0 += 1
            return rxns_idx0
        end

        return true
    end

    ## ..- -- .- . --
    # cobra fields
    M = length(mets_idxmap)
    Ne = length(mets_e) # exchange reactions
    N = length(rxns_idxmap) 
    @show M, N

    mets = Vector{String}(undef, M)
    rxns = Vector{String}(undef, N)
    subSystems = Vector{String}(undef, N)
    rxnNames = Vector{String}(undef, N)
    lb = Vector{Float64}(undef, N)
    ub = Vector{Float64}(undef, N)
    S = spzeros(Float64, M, N) #TODO/TAI make sparse
    b = zeros(M)
    c = zeros(N)

    # ast reactions
    for (rxn, ast) in pannet_ast
        isempty(rxn) && continue
        rxni = rxns_idxmap[rxn]

        rxns[rxni] = rxn
        subSystems[rxni] = ast["subSys"]
        rxnNames[rxni] = ismissing(ast["rxnName"]) ? "" : ast["rxnName"]
        lb[rxni] = ast["lb"]
        ub[rxni] = ast["ub"]

        for (met, s) in ast["stoi"]
            isempty(met) && continue
            meti = mets_idxmap[met]
            mets[meti] = met
            S[meti, rxni] = s
        end
    end

    # rxn[344]: EX_glc__D_e (D-Glucose exchange)
    # subsys: nothing
    # lb: -10.0, ub: 0.0
    # (-1.0) glc__D_e <== 
    for met_e in mets_e
        meti = mets_idxmap[met_e]
        rxn = string("EX_", met_e)
        rxni = rxns_idxmap[rxn]

        rxns[rxni] = rxn
        subSystems[rxni] = "Exchange"
        rxnNames[rxni] = string(met_e, " exchange reaction")
        lb[rxni] = 0
        ub[rxni] = 1000
        S[meti, rxni] = -1.0
    end

    #MARK: MetNet
    global net0 = MetNet(;
        S, b, c, lb, ub, mets, rxns, 
        subSystems, rxnNames
    )

    # copy iML1515 bounds
    global iML1515_net0 = pull_net("iML1515")
    for rxn in reactions(iML1515_net0)
        if hasrxnid(net0, rxn)
            il, iu = bounds(iML1515_net0, rxn)
            bounds!(net0, rxn, il, iu)
        else
            # summary(iML1515_net0, rxn)
        end
    end

    # MARK: Medium
    # Close sinks
    for rxn in net0.rxns
        startswith(rxn, "sink_") || continue
        bounds!(net0, rxn, 0.0, 0.0)
    end

    # Maintenance
    # TODO: Search value in literature
    bounds!(net0, "ATPM", 10.0, 100000.0)
    # bounds!(net0, "ATPM", 610.0, 100000.0)
    
    
    # Close all intake/ Open outtake
    for rxn in net0.rxns
        startswith(rxn, "EX_") || continue
        bounds!(net0, rxn, 0.0, 1000.0)
        # bounds!(net0, rxn, -1000.0, 1000.0)
    end
    
    bounds!(net0,"EX_glc__D_e", -10.0, 0.0)
    # bounds!(net0,"EX_o2s_e", -1000.0, 1000.0)
    bounds!(net0,"EX_o2_e", -1000.0, 1000.0)

    # bounds!(net0,"NDPK1", 0.0, 0.0)
    # bounds!(net0,"BGLA", 0.0, 0.0)
    # bounds!(net0,"MANE", 0.0, 0.0)
    # bounds!(net0,"GALM1", 0.0, 0.0)
    # bounds!(net0,"AMALT4", 0.0, 0.0)

    # bounds!(net0,"SUCCtex", 0.0, 0.0)
    # bounds!(net0,"SUCFUMtpp", 0.0, 0.0)
    # bounds!(net0,"TARTRt7pp", 0.0, 0.0)
    # bounds!(net0,"SUCTARTtpp", 0.0, 0.0)
    # bounds!(net0,"O2tex", -1000.0, 100.0)

    # bounds!(net0,"EX_h2o2_e", 0.0, 1000.0)
    # bounds!(net0,"EX_glcr_e", 0.0, 1000.0)
    # bounds!(net0, "EX_metglcur_e", 0.0, 1000.0)
    # bounds!(net0,"EX_glcur_e", 0.0, 1000.0)
    # bounds!(net0,"EX_2dhglcn_e", 0.0, 1000.0)
    # bounds!(net0,"EX_udpglcur_e", 0.0, 1000.0)
    # bounds!(net0,"EX_25dkglcn_e", 0.0, 1000.0)
    # bounds!(net0,"EX_udpglcur_e", 0.0, 1000.0)
    # bounds!(net0,"EX_2ddglcn_e", 0.0, 1000.0)
    # bounds!(net0,"EX_fruur_e", 0.0, 1000.0)
    # bounds!(net0,"EX_glcn_e", 0.0, 1000.0)
    # bounds!(net0,"EX_glcur1p_e", 0.0, 1000.0)

    bounds!(net0,"EX_mg2_e", -1000.0, 1000.0)
    bounds!(net0,"EX_so3_e", -1000.0, 1000.0)
    bounds!(net0,"EX_cl_e", -1000.0, 1000.0)
    bounds!(net0,"EX_ca2_e", -1000.0, 1000.0)
    bounds!(net0,"EX_cu_e", -1000.0, 1000.0)
    bounds!(net0,"EX_k_e", -1000.0, 1000.0)
    bounds!(net0,"EX_zn2_e", -1000.0, 1000.0)
    # bounds!(net0,"EX_h2_e", -1000.0, 1000.0)
    # bounds!(net0, "EX_co2_e", -1000.0, 1000.0)
    bounds!(net0, "EX_fe2_e", -1000.0, 1000.0)
    bounds!(net0, "EX_nh4_e", -1000.0, 1000.0)
    bounds!(net0, "EX_mn2_e", -1000.0, 1000.0)
    bounds!(net0, "EX_pi_e", -1000.0, 1000.0)
    bounds!(net0, "EX_cobalt2_e", -1000.0, 1000.0)

    bounds!(net0, "EX_h2o_e", -1000.0, 1000.0)
    bounds!(net0, "EX_h_e", -1000.0, 1000.0)
    
    # bounds!(net0,"EX_so4_e", -1000.0, 1000.0)
    # bounds!(net0, "EX_na1_e", -1000.0, 1000.0)
    # bounds!(net0, "EX_k_e", -1000.0, 1000.0)
    # bounds!(net0, "EX_cu2_e", -1000.0, 1000.0)
    # bounds!(net0,"EX_n2o_e", -1000.0, 1000.0)
    
    # extras
    extras!(net0, "BIOM", "Growth")

    extras!(net0, "EX_GLC", "EX_glc__D_e")
    extras!(net0, "EX_NH4", "EX_nh4_e")
    extras!(net0, "EX_GLU", "EX_glu__L_e")
    extras!(net0, "EX_O2", "EX_o2_e")
    extras!(net0, "EX_CO2", "EX_co2_e")
    # extras!(net, "ATPM", "ATPM")

    linear_weights!(net0, "Growth", 1.0)

    #MARK: FBA test
    # opm = FBAOpModel(net0, Ipopt.Optimizer)
    global opm = FBAOpModel(net0, Clp.Optimizer)
    
    set_linear_obj!(opm, "Growth", MAX_SENSE)
    optimize!(opm)
    biom = solution(opm, "Growth")
    @show biom


    global open_rxns = []
    fn = joinpath(@__DIR__, "fba.md")
    open(fn, "w") do io
        println(io, "-"^20)
        println(io, "EX_reactions")
        println(io, "\n\n\n")

        for rxn in net0.rxns
            startswith(rxn, "EX_") || continue
            sol = solution(opm, rxn)
            abs(sol) > 1e2 || continue
            MetX.MetXGEMs._print_rxn_summary(io, net0, rxn)
            println(io, "sol: ", sol)
            println(io)
        end
        
        println(io, "\n\n\n")
        println(io, "-"^20)
        println(io, "Inner_reactions")
        for rxn in net0.rxns
            startswith(rxn, "EX_") && continue
            sol = solution(opm, rxn)
            abs(sol) > 1e2 || continue
            MetX.MetXGEMs._print_rxn_summary(io, net0, rxn)
            println(io, "sol: ", sol)
            println(io)
        end
    end
    return 

    # bounds!(opm, "Growth", biom - 1e-3, biom)
    
    # set_linear_obj!(opm, "EX_glc__D_e", MIN_SENSE)
    # optimize!(opm)
    # biom = solution(opm, "Growth")
    # ex_glc = solution(opm, "EX_glc__D_e")
    # @show biom
    # @show ex_glc
    # bounds!(opm, "EX_glc__D_e", ex_glc, ex_glc * 0.5)
    
    # set_v2_obj!(opm, MAX_SENSE)
    # optimize!(opm)
    # @show solution(opm, "Growth")

    global open_rxns = []
    for rxn in net0.rxns
        startswith(rxn, "EX_") || continue
        sol = solution(opm, rxn)
        sol < 0 || continue
        @show rxn, sol
        push!(open_rxns, rxn)
    end
    
    # net0

end

## ..- -- .- . -- -- - . .. .. -. -.- 
let
    global opm = FBAOpModel(net0, Clp.Optimizer)
   
    sol0 = nothing
    _diffs = []
    for atpm in 1:20:100
        bounds!(opm, "ATPM", atpm, 100000.0)

        set_linear_obj!(opm, "Growth", MAX_SENSE)
        optimize!(opm)
        biom = solution(opm, "Growth")
        
        sol = solution(opm)
        if isnothing(sol0)
            sol0 = sol
            continue
        end
        
        if isempty(_diffs)
            _diffs = sol0 - sol
            continue
        end

        _diffs = [_diffs sol0 - sol]


        @show biom
    end

    f = Figure()
    ax = Axis(f[1,1])
    for (ri, r) in enumerate(eachrow(_diffs))
        mean(abs, r) > 1.3e3 || continue
        lines!(ax, eachindex(r), r)
        
        summary(net0, reactions(net0, ri))
    end
    f

end


## ..- -- .- . -- -- - . .. .. -. -.- 
let
    f = Figure()
    ax = Axis(f[1,1])
    M = rand(10,5)
    for r in eachrow(M)
        lines!(ax, eachindex(r), r)
    end
    f
end


## ..- -- .- . -- -- - . .. .. -. -.- 
let
    for rxn in open_rxns
        println("bounds!(net0,", repr(rxn), ", -1000.0, 1000.0)")
    end
end

## ..- -- .- . -- -- - . .. .. -. -.- 
let
    global iJR904_net0 = pull_net("iJR904")
    for rxn in iJR904_net0.rxns
        i9_l, i9_u = bounds(iJR904_net0, rxn)
        l, u = bounds(net0, rxn)
        i9_l != l || continue
        summary(iJR904_net0, rxn)
        summary(net0, rxn)
        break
    end
end

## ..- -- .- . -- -- - . .. .. -. -.- 
let
    global iJR904_net0 = pull_net("iJR904")
    bounds!(iJR904_net0, "EX_o2_e", -0.0, 1000.0)

    #MARK: FBA test
    opm = FBAOpModel(iJR904_net0, Clp.Optimizer)
    optimize!(opm)
    @show solution(opm, extras(iJR904_net0, "BIOM"))
    opm

    global open_rxns = []
    for rxn in iJR904_net0.rxns
        startswith(rxn, "EX_") || continue
        sol = solution(opm, rxn)
        sol < 0 || continue
        @show rxn, sol
        push!(open_rxns, rxn)
    end
end

## ..- -- .- . -- -- - . .. .. -. -.- 
#MARK: replicate fig 3a
let
    all_rxns = rxninfo_tsv[:, "Reaction ID"]
    all_rxn_names = rxninfo_tsv[:, "Reaction Name"]
    idx = rand(1:3342)
    @show idx
    rxninfo_tsv[:, "Reaction"][idx]

    rxn_ast["rxn"] = rxn
    # @show rxn
    rxnName = rxninfo_tsv[:, "Reaction Name"][ridx]
    rxn_ast["rxnName"] = rxnName
    # @show rxn_name
    rxnEq = rxninfo_tsv[:, "Reaction"][ridx]

    JSON.print(rxn_ast, 4)
end