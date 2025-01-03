@time begin
    using JSON
    using CSV
    using DataFrames
end

let
    csvdir = joinpath(@__DIR__, "proc")
    
    # csvfile = joinpath(csvdir, "Batch-elife-79815-supp1-v2.csv")
    # csvfile = joinpath(csvdir, "Chemostat-elife-79815-supp1-v2.csv")
    csvfile = joinpath(csvdir, "Batch-elife-79815-supp2-v2.csv")
    
    
    df = CSV.read(csvfile, DataFrame; delim = ',', stringtype=String)
    
    dict = Dict()
    for colname in names(df)
        dict[colname] = Dict(
            "data" => df[:, colname]
        )
    end
    
    jsonfile = replace(csvfile, ".csv" => ".json")
    @show jsonfile
    open(jsonfile, "w") do io
        println(io, JSON.json(dict))
    end

end