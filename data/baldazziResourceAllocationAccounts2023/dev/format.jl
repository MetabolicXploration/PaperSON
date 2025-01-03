@time begin
    using DataFrames
    using CSV
end

## --.. - . -.. -.- .-.. - .- -.-- -- . . . . .- - 
# fill column
function _down_complete!(fun::Function, df::DataFrame, col)
    _last = missing
    for row in eachrow(df)
        if fun(row[col]) 
            row[col] = _last
        else
            _last = row[col]
        end
    end
    return df
end
_down_complete!(df::DataFrame, col) = _down_complete!(ismissing, df, col)

## --.. - . -.. -.- .-.. - .- -.-- -- . . . . .- - 
# supp1 and supp2
# 
let
    # read
    for file in [
            "Chemostat-elife-79815-supp1-v2.tsv",
            "Batch-elife-79815-supp1-v2.tsv",
            "Batch-elife-79815-supp2-v2.tsv"
        ]

        @show file
        
        rawfile = joinpath(@__DIR__, "raw", file)
        df = CSV.read(rawfile, DataFrame; delim = '\t', stringtype=String)
        
        # format col names
        # make homogeneous names
        for colname0 in names(df)
            colname = strip(colname0)
            colname = replace(colname, 
                "Flux unit" => "Rate unit",
                "Biomass Yield (gDW/g_glucose)" => "Biomass Yield (gDW/g_glc)"
            )
            rename!(df, colname0 => colname)
        end
        
        # complete
        _down_complete!(df, "Rate unit")
        _down_complete!(df, "Reference")
        _down_complete!(df, "PMID")

        # make homogeneous values
        for row in eachrow(df)
            row["Rate unit"] = replace(row["Rate unit"], " "=>"")
        end
        
        # write
        procfile = joinpath(@__DIR__, "proc", basename(rawfile))
        mkpath(dirname(procfile))
        CSV.write(procfile, df)
        
        # test load
        df = CSV.read(procfile, DataFrame; delim = ',')
        display(first(df, 5))
    end
end
## --.. - . -.. -.- .-.. - .- -.-- -- . . . . .- - 