## --.. - . -.. -.- .-.. - .- -.-- -- . . . . .- - 
@time begin
    using DataFrames
    using CSV
    using JSON
end

## --.. - . -.. -.- .-.. - .- -.-- -- . . . . .- - 
let
    tsv_file = joinpath(@__DIR__, "rstb20210236_si_002-GPRs.tsv")

    global df = CSV.read(tsv_file, DataFrame)

    # coalesce.(df[:, 10], 0)
    global dict = Dict{String, Vector}()
    
    _fill = Dict(
        "Reaction ID" => "", 
        "Reaction" => "",
        "Reaction Name"  => "",
        "Reaction Count" => 0,
        "Core/Accessory" => "",
        "Subsystem" => "",
        "default" => ""
    )
    
    for colid in names(df)
        # @show colid
        vals = df[:, colid]
        val0 = get(_fill, colid, _fill["default"])
        dict[colid] = coalesce.(df[:, colid], val0)
        # @show val0
    end
    
    json_file = joinpath(@__DIR__, "..", "raw.rstb20210236_si_002-GPRs.json_")
    open(json_file, "w") do io
        print(io, "{")
        print(io, "\n")
        for id in names(df)
            dat = dict[id]
            print(io, "    ")
            JSON.print(io, id)
            print(io, " : ")
            JSON.print(io, dat)
            print(io, ",")
            print(io, "\n\n")
        end
        print(io, "\n")
        print(io, "}")
    end
end

## --.. - . -.. -.- .-.. - .- -.-- -- . . . . .- - 
let
    json_file = joinpath(@__DIR__, "..", "raw.rstb20210236_si_002-GPRs.json_")
    dict = JSON.parsefile(json_file)
end