@time begin
    using CSV
    using DataFrames
    using JSON
    using CairoMakie
end

## - --- -- . .- -.-. -.- .-.-. 
# TODO/ Add to folder README

## - --- -- . .- -.-. -.- .-.-. 
# MARK: load
_tonumeric(vec::Vector, fill = NaN) = [isa(x, Number) ? x : fill for x in vec]

let
    json_file = joinpath(@__DIR__, "raw.json")
    global rawData = JSON.parsefile(json_file)
    json_file = joinpath(@__DIR__, "raw.msystems.00625-20-sd001.xlsx.json")
    global sd001Data = JSON.parsefile(json_file)

    global sheet4Data = sd001Data["data"]["Sheet 4"]["data"]
    global fig3Data = rawData["data"]["Fig 3"]
    nothing
end

## - --- -- . .- -.-. -.- .-.-. 
# MARK: Fig 3

## - --- -- . .- -.-. -.- .-.-. 
# A
let 
    f = Figure()
    ax = Axis(f[1,1];
        title = "Biomass Production",
        xlabel = fig3Data["Biomass Production"]["xlabel"],
        ylabel = fig3Data["Biomass Production"]["ylabel"],
        limits = (0.0, 10.0, -0.05, 0.8)
    )
    
    xs = _tonumeric(sheet4Data["Glucose uptake rate"]["vals"], NaN)
    ys = _tonumeric(sheet4Data["Growth rate"]["vals"], NaN)
    scatter!(ax, xs, ys, label = "experimental")
    
    
    xs = _tonumeric(fig3Data["Biomass Production"]["PAM Simulations"]["x"])
    ys = _tonumeric(fig3Data["Biomass Production"]["PAM Simulations"]["y"])
    lines!(ax, xs, ys)
    f
end

## - --- -- . .- -.-. -.- .-.-. 
# B
let 
    f = Figure()
    ax = Axis(f[1,1];
        title = "Acetate Secretion",
        xlabel = fig3Data["Acetate Secretion"]["xlabel"],
        ylabel = fig3Data["Acetate Secretion"]["ylabel"],
        limits = (0.0, 10.0, -0.5, 8.0)
    )
    
    xs = _tonumeric(sheet4Data["Glucose uptake rate"]["vals"], NaN)
    ys = _tonumeric(sheet4Data["Acetate secretion rate"]["vals"], NaN)
    scatter!(ax, xs, ys, label = "experimental")
    
    
    xs = _tonumeric(fig3Data["Acetate Secretion"]["PAM Simulations"]["x"])
    ys = _tonumeric(fig3Data["Acetate Secretion"]["PAM Simulations"]["y"])
    lines!(ax, xs, ys)
    f
end

## - --- -- . .- -.-. -.- .-.-. 
# C
let 
    f = Figure()
    ax = Axis(f[1,1];
        title = "O2 Production",
        xlabel = fig3Data["O2 Production"]["xlabel"],
        ylabel = fig3Data["O2 Production"]["ylabel"],
        limits = (0.0, 10.0, -0.5, 20.0)
    )
    
    xs = _tonumeric(sheet4Data["Glucose uptake rate"]["vals"], NaN)
    ys = _tonumeric(sheet4Data["O2 uptake rate"]["vals"], NaN)
    scatter!(ax, xs, ys, label = "experimental")
    
    
    xs = _tonumeric(fig3Data["O2 Production"]["PAM Simulations"]["x"])
    idx = sortperm(xs)
    ys = _tonumeric(fig3Data["O2 Production"]["PAM Simulations"]["y"])
    # clipboard(xs[idx])
    # clipboard(ys[idx])
    # lines!(ax, xs[idx], ys[idx])
    lines!(ax, xs, ys)
    f
end

## - --- -- . .- -.-. -.- .-.-. 
# D
let 
    f = Figure()
    ax = Axis(f[1,1];
        title = "CO2 Production",
        xlabel = fig3Data["CO2 Production"]["xlabel"],
        ylabel = fig3Data["CO2 Production"]["ylabel"],
        limits = (0.0, 10.0, -0.5, 20.0)
    )
    
    xs = _tonumeric(sheet4Data["Glucose uptake rate"]["vals"], NaN)
    ys = _tonumeric(sheet4Data["CO2 production rate"]["vals"], NaN)
    scatter!(ax, xs, ys, label = "experimental")
    
    xs = _tonumeric(fig3Data["CO2 Production"]["PAM Simulations"]["x"])
    idx = sortperm(xs)
    ys = _tonumeric(fig3Data["CO2 Production"]["PAM Simulations"]["y"])
    lines!(ax, xs, ys)
    f
end

## - --- -- . .- -.-. -.- .-.-. 