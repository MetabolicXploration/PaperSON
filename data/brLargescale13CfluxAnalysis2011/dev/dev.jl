@time begin
    using CSV
    using DataFrames
    using JSON
    using OrderedCollections
end

function _foreach_key(f::Function, root_dict, key0)
    for (key, dat) in root_dict
        # base
        if key == key0
            f(key, root_dict) === :break && return :break
        end

        # recursive
        if dat isa AbstractDict
            _foreach_key(f, dat, key0) === :break && return :break
        end
    end
end

let
    tsv_file = joinpath(@__DIR__, "dev.tsv")
    df = CSV.read(tsv_file, DataFrame; header = 0)

    json_mask = joinpath(@__DIR__, "mask.json")
    mask = JSON.parsefile(json_mask; dicttype = OrderedDict)

    R, C = size(df)
    @show (R, C)
    ri = 1
    ci = 1
    # (ri, ci) = (93, 37)
    # It will go row by row
    _foreach_key(mask, "val") do k, dict
        
        @show (ri, ci)
        dict["val"] = df[ri, ci]
        ci += 1
        if ci > C
            ci = 1
            ri += 1
            # ri == 3 && return :break
        end
        
        return :continue
    end
    @show (R, C)
    
    
    # return
    
    json_out = joinpath(@__DIR__, "out.json")
    open(json_out, "w") do io
        JSON.print(io, mask, 4)
    end
end