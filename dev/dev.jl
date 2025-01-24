@time begin
    using JSON
    using CSV
    using DataFrames
    using CairoMakie
end

## ..- - - - --- ..  . ... -. -. -.- -.. .- - ... -
# (-0.05) 5mthf_c + (-5.0e-5) accoa_c + (-0.488) ala__L_c + (-0.001) amp_c + (-0.281) arg__L_c + (-0.229) asn__L_c + (-0.229) asp__L_c + (-45.7318) atp_c + (-0.000129) clpn_EC_c + (-6.0e-6) coa_c + (-0.126) ctp_c + (-0.087) cys__L_c + (-0.0247) datp_c + (-0.0254) dctp_c + (-0.0254) dgtp_c + (-0.0247) dttp_c + (-1.0e-5) fad_c + (-0.25) gln__L_c + (-0.25) glu__L_c + (-0.582) gly_c + (-0.154) glycogen_c + (-0.203) gtp_c + (-45.5608) h2o_c + (-0.09) his__L_c + (-0.276) ile__L_c + (-0.428) leu__L_c + (-0.0084) lps_EC_c + (-0.326) lys__L_c + (-0.146) met__L_c + (-0.00215) nad_c + (-5.0e-5) nadh_c + (-0.00013) nadp_c + (-0.0004) nadph_c + (-0.001935) pe_EC_c + (-0.0276) peptido_EC_c + (-0.000464) pg_EC_c + (-0.176) phe__L_c + (-0.21) pro__L_c + (-5.2e-5) ps_EC_c + (-0.035) ptrc_c + (-0.205) ser__L_c + (-0.007) spmd_c + (-3.0e-6) succoa_c + (-0.241) thr__L_c + (-0.054) trp__L_c + (-0.131) tyr__L_c + (-0.003) udpg_c + (-0.136) utp_c + (-0.402) val__L_c ==> (45.5608) adp_c + (45.56035) h_c + (45.5628) pi_c + (0.7302) ppi_c

## ..- - - - --- ..  . ... -. -. -.- -.. .- - ... -
let
    datdir = "/Users/Pereiro/University/Research/CODE/PaperSON/data/baldazziResourceAllocationAccounts2023"

    jsonfile = joinpath(datdir, "raw.elife-79815-supp1-v2.json")
    global d = JSON.parsefile(jsonfile)

    # filter rows
    idx0 = findall(d["data"]["Batch"]["Rate unit"]["data"]) do el
        el == "mmol/gDW/h"
    end

    id = "Glucose uptake rate"
    # id = "Growth rate (1/h)"
    # id = "Biomass Yield (gDW/g_glc)"
    # id = "Yield"
    dat = d["data"]["Batch"][id]["data"]
    idx1 = findall(!isnothing, dat)
    
    idx = intersect(idx0, idx1)

    f = Figure()
    ax = Axis(f[1,1]; xlabel = id, ylabel = "count")
    hist!(ax, dat[idx]; bins = 20)
    f
end

## ..- - - - --- ..  . ... -. -. -.- -.. .- - ... -
let
    datdir = "/Users/Pereiro/University/Research/CODE/PaperSON/data/baldazziResourceAllocationAccounts2023"

    jsonfile = joinpath(datdir, "raw.elife-79815-supp1-v2.json")
    global d = JSON.parsefile(jsonfile)

    id1 = "Growth rate (1/h)"
    id2 = "Glucose uptake rate"

    f = Figure()
    ax1 = Axis(f[1,1]; 
        xlabel = id1, ylabel = id2,
        limits = (0.0, nothing, 0.0, nothing)
    )
    ax2 = Axis(f[2,1];
        xlabel = "Yield (gWD/Cmmol)", 
        ylabel = "count",
        limits = (0.0, 0.03, 0.0, nothing)
    )

    for (bins, color, cultype) in [
            # (25, :blue, "Batch"), 
            (10, :red, "Chemostat"), 
        ]
        # filter rows
        idx0 = findall(d["data"][cultype]["Rate unit"]["data"]) do el
            el == "mmol/gDW/h"
        end

        # id1 = "Glucose uptake rate"
        # id1 = "Biomass Yield (gDW/g_glc)"
        dat1 = d["data"][cultype][id1]["data"]
        @show length(dat1)
        idx1 = findall(!isnothing, dat1)

        # id2 = "Biomass Yield (gDW/g_glc)"
        # Carbon content per gDW E. coli
        dat2 = d["data"][cultype][id2]["data"]
        @show length(dat2)
        idx2 = findall(!isnothing, dat2)

        idx = intersect(idx0, idx1, idx2)

        
        scatter!(ax1, dat1[idx], dat2[idx]; color)
        Yg = dat1[idx] ./ (dat2[idx] .* 6)
        @show extrema(Yg)
        hist!(ax2, Yg;  
            color, bins, 
            normalization = :probability, 
            # alpha = 0.2
        )
    end
    f
end

## ..- - - - --- ..  . ... -. -. -.- -.. .- - ... -
let
    datdir = "/Users/Pereiro/University/Research/CODE/PaperSON/data/baldazziResourceAllocationAccounts2023"

    jsonfile = joinpath(datdir, "raw.elife-79815-supp2-v2.json")
    global d = JSON.parsefile(jsonfile)

    id1 = "Growth rate (1/h)"
    id2 = "Glycerol uptake rate"

    f = Figure()
    ax1 = Axis(f[1,1]; 
        xlabel = id1, ylabel = id2,
        limits = (0.0, nothing, 0.0, nothing)
    )
    ax2 = Axis(f[2,1];
        xlabel = "Yield (gWD/Cmmol)", 
        ylabel = "count",
        limits = (0.0, 0.03, 0.0, nothing)
    )

    # filter rows
    idx0 = findall(d["data"]["Rate unit"]["data"]) do el
        el == "mmol/gDW/h"
    end

    # id1 = "Glucose uptake rate"
    # id1 = "Biomass Yield (gDW/g_glc)"
    dat1 = d["data"][id1]["data"]
    @show length(dat1)
    idx1 = findall(!isnothing, dat1)

    # id2 = "Biomass Yield (gDW/g_glc)"
    # Carbon content per gDW E. coli
    dat2 = d["data"][id2]["data"]
    @show length(dat2)
    idx2 = findall(!isnothing, dat2)

    idx = intersect(idx0, idx1, idx2)

    
    scatter!(ax1, dat1[idx], dat2[idx])
    hist!(ax2, dat1[idx] ./ (dat2[idx] .* 3); 
        bins = 10
    )
    f
end

## ..- - - - --- ..  . ... -. -. -.- -.. .- - ... -
let
    datdir = "/Users/Pereiro/University/Research/CODE/PaperSON/data/baldazziResourceAllocationAccounts2023"

    jsonfile = joinpath(datdir, "raw.elife-79815-supp1-v2.json")
    global d = JSON.parsefile(jsonfile)

    # filter rows
    idx0 = findall(d["data"]["Batch"]["Rate unit"]["data"]) do el
        el == "mmol/gDW/h"
    end

    # id1 = "Glucose uptake rate"
    # id1 = "Biomass Yield (gDW/g_glc)"
    id1 = "Growth rate (1/h)"
    dat1 = d["data"]["Batch"][id1]["data"]
    @show length(dat1)
    idx1 = findall(!isnothing, dat1)

    # id2 = "Biomass Yield (gDW/g_glc)"
    id2 = "Yield"
    dat2 = d["data"]["Batch"][id2]["data"]
    @show length(dat2)
    idx2 = findall(!isnothing, dat2)

    idx = intersect(idx0, idx1, idx2)

    f = Figure()
    ax = Axis(f[1,1]; 
        xlabel = id1, ylabel = id2,
        limits = (0.0, nothing, 0.0, nothing)
    )
    scatter!(ax, dat1[idx], dat2[idx])
    f
end

## ..- - - - --- ..  . ... -. -. -.- -.. .- - ... -
let
    datdir = "/Users/Pereiro/University/Research/CODE/PaperSON/data/baldazziResourceAllocationAccounts2023"

    jsonfile = joinpath(datdir, "raw.elife-79815-supp1-v2.json")
    global d = JSON.parsefile(jsonfile)
    chtat = d["data"]["Chemostat"]["Biomass Yield (gDW/g_glc)"]["data"]
    chtat = filter(!isnothing, chtat)
    # chtat ./= 6
    batch = d["data"]["Batch"]["Biomass Yield (gDW/g_glc)"]["data"]
    batch = filter(!isnothing, batch)
    # batch ./= 6

    f = Figure()
    ax = Axis(f[1,1]; xlabel = "yield", ylabel = "count")
    # hist!(ax, chtat; color = :blue)
    # hist!(ax, batch; color = :red)
    hist!(ax, [batch; chtat])
    f
end

## ..- - - - --- ..  . ... -. -. -.- -.. .- - ... -
let
    datdir = "/Users/Pereiro/University/Research/CODE/PaperSON/data/baldazziResourceAllocationAccounts2023"

    jsonfile = joinpath(datdir, "raw.elife-79815-supp2-v2.json")
    global d = JSON.parsefile(jsonfile)
    chtat = d["data"]["Yield"]["data"]
    chtat = filter(!isnothing, chtat)
    # chtat ./= 6

    f = Figure()
    ax = Axis(f[1,1]; xlabel = "yield", ylabel = "count")
    hist!(ax, chtat; color = :blue)
    # hist!(ax, batch; color = :red)
    # hist!(ax, [batch; chtat])
    f
end

## ..- - - - --- ..  . ... -. -. -.- -.. .- - ... -