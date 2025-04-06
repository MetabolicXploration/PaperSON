@time begin
    using XLSX
end

## .- - .- . -.- -... .- -- - - ---- . . .
# TODO: create a function that for a given allow you to extract sections of a xlsx file to an tsv file. You can creates a mapping interface (sheet, xlsx_range) -> tsv_range. 
let
    sep = "\t"
    src_path = "/Users/Pereiro/University/Research/CODE/PaperSON/data/alterProteomeRegulationPatterns2021/dev/msystems.00625-20-sd001.xlsx"
    
    xf = XLSX.readxlsx(src_path)

    range_pool = Dict{String, Any}()
    # range_pool["Sheet 1"] = "A6:F2849"
    
    for (sheetname) in XLSX.sheetnames(xf)
        @show sheetname
        out_path = string(src_path, "-", sheetname, ".tsv")
        @show out_path
        open(out_path, "w") do io
            sh = xf[sheetname]
            range = get(range_pool, sheetname, Colon())
            mat = sh[range]
            @show size(mat)
            for row in eachrow(mat)
                row_str = join(string.(row), "\t")
                println(io, row_str)
            end
            println(io)
        end
    end
end


## .- - .- . -.- -... .- -- - - ---- . . .