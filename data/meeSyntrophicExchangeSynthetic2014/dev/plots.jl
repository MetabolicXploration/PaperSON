@time begin
    using CairoMakie
    using CSV
    using JSON
    using DataFrames
end

## ..- .. - - -. .. .- - - ..- -. -..- -.- .
PAPERSON_DIR = "/Users/pereiro/University/CODE/PaperSON"
PAPERSON_DATA_DIR = joinpath(PAPERSON_DIR, "data")

## ..- .. - - -. .. .- - - ..- -. -..- -.- .
let
    MEE_DIR = joinpath(PAPERSON_DATA_DIR, 
        "meeSyntrophicExchangeSynthetic2014"
    )
    
    raw_path = joinpath(MEE_DIR, "raw.json")
   
    raw = JSON.parsefile(raw_path)
    scv_name = raw["data"]["pnas.1405641111.sd03"]["2-member"]["data"]["filename"]
    scv_file = joinpath(MEE_DIR, scv_name)

    global df = CSV.read(scv_file, DataFrame; 
        header=2
    )

    names(df)
end


## ..- .. - - -. .. .- - - ..- -. -..- -.- .
let
    # df[!, "Co-culture ID"]
    # df[!, "Strain 1"]
    # df[!, "Strain 1\nfold growth\n(T84)"]
    # df[!, "coop coeff\nc_12"]
    # df[!, "Strain 2"]
    # df[!, "Strain 2\nfold growth\n(T84)"]
    # df[!, "coop coeff\nc_21"]

    # st_ids = unique!(df[:, "Strain 1"])
    global DAT = Dict()
    for row in eachrow(df)
        si1 = String(row["Strain 1"])
        si2 = String(row["Strain 2"])
        f1 = row["Strain 1\nfold growth\n(T84)"]
        f2 = row["Strain 2\nfold growth\n(T84)"]
        dat1 = get!(DAT, si1, Float64[])

        x0 = 1e7 / 2 # cells / ml
        dx1 = x0 * f1 - x0
        dx2 = x0 * f2 - x0
        
        # Syntrophy only
        dx1 > 0 || continue
        dx2 > 0 || continue

        dat1 = get!(DAT, si1, Float64[])
        push!(dat1, log10(dx1) ./ log10(dx2))
        dat2 = get!(DAT, si2, Float64[])
        push!(dat2, log10(dx2) ./ log10(dx1))
    end

end

## ..- .. - - -. .. .- - - ..- -. -..- -.- .
let
    # Rx7CRkXbRFbJdC3b
    f = Figure()
    ax = Axis(f[1,1]; 
        title = "Syntrophic cultures",
        ylabel = "trade yield (log/log)",
        xlabel = "parner index (unsorted)"
    )
    for (si1, yields) in DAT
        
        ys = sort(yields)
        xs = eachindex(ys)
        # barplot!(ax, xs, ys; 
        #     label = si1, 
        #     alpha = 0.4
        # )
        scatter!(ax, xs, ys; 
            markersize = 16
        )
        lines!(ax, xs, ys)
    end
    f
end

## ..- .. - - -. .. .- - - ..- -. -..- -.- .
## ..- .. - - -. .. .- - - ..- -. -..- -.- .
let    


    s1 = df[!, "Strain 1\nfold growth\n(T84)"]
    s2 = df[!, "Strain 2\nfold growth\n(T84)"]

    f = Figure()
    ax = Axis(f[1,1])

    ys = s1 ./ s2
    sort!(ys)
    xs = eachindex(ys)
    scatter!(ax, xs, ys)
    f
end

## ..- .. - - -. .. .- - - ..- -. -..- -.- .
## ..- .. - - -. .. .- - - ..- -. -..- -.- .
# df = CSV.read(raw_path, DataFrame; 
#     validate = true, 
#     strict = false,
#     skipto = 5,
# )
## ..- .. - - -. .. .- - - ..- -. -..- -.- .
## ..- .. - - -. .. .- - - ..- -. -..- -.- .