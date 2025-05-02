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
    global fig1Data = rawData["data"]["Fig 1"]
    nothing
end


## - --- -- . .- -.-. -.- .-.-. 
# MARK: Figure 1

function _plot_fig1!(ax, key)
    dat = fig1Data["data"][key]
    ax.xlabel = dat["x.label"]
    ax.ylabel = dat["y.label"]
    scatter!(ax, dat["x.vals"], dat["y.vals"])
end

## - --- -- . .- -.-. -.- .-.-. 
let
    f = Figure()
    ax = Axis(f[1:3,1:3])
    _plot_fig1!(ax, "Experimental.Glucose")
    f
end

## - --- -- . .- -.-. -.- .-.-. 
let
    f = Figure()
    ax = Axis(f[1:3,1:3])
    _plot_fig1!(ax, "Experimental.Acetate")
    f
end

## - --- -- . .- -.-. -.- .-.-. 
let
    f = Figure()
    ax = Axis(f[1:3,1:3])
    _plot_fig1!(ax, "Experimental.Biomass")
    f
end