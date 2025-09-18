#!/usr/bin/env julia
# xlsx2csv.jl
#
# Usage examples:
#   julia xlsx2csv.jl input.xlsx                 # first sheet → input.csv
#   julia xlsx2csv.jl input.xlsx --sheet=Data    # named sheet → input_Data.csv
#   julia xlsx2csv.jl input.xlsx --sheet=2       # 2nd sheet → input_Sheet2.csv
#   julia xlsx2csv.jl input.xlsx --range=A1:D20  # rectangular range on first sheet
#   julia xlsx2csv.jl input.xlsx --all           # all sheets → one CSV per sheet
#   julia xlsx2csv.jl input.xlsx --out=out.csv   # explicit output path (single sheet only)

## ---.... .. - . .-  . -. . . - .  . . - -
using XLSX
using DataFrames
using CSV

## ---.... .. - . .-  . -. . . - .  . . - -
# MARK: Utils
function sheet_to_matrix(file::AbstractString; sheet=1, range=nothing)
    XLSX.openxlsx(file) do xf
        sname = sheet isa Integer ? XLSX.sheetnames(xf)[sheet] : String(sheet)
        sh = xf[sname]

        # This is the right way to get the sheet's used range
        used_ref = XLSX.used_range(sh)

        if range === nothing
            return XLSX.readdata(sh, used_ref)
        else
            return XLSX.readdata(sh, range)
        end
    end
end


## ---.... .. - . .-  . -. . . - .  . . - -
let
    data_folder = "/Users/pereiro/University/CODE/PaperSON/data/meeSyntrophicExchangeSynthetic2014/dev"
    
    xlsx_name = "pnas.1405641111.sd03.xlsx"
    xlsx_file = joinpath(data_folder, xlsx_name)
    

    # Open the Excel file
    xlsx = XLSX.readxlsx(xlsx_file)

    # Select a specific sheet by name or index
    # sheet_name = "2-member"
    # sheet_range = "A:G"

    sheet_name = "3-member"
    sheet_range = "A:H"

    # sheet_name = "14&13-member mean"
    # sheet_range = "A:AO"
    sheet = xlsx[sheet_name]  # or xlsx[1] for the first sheet

    csv_file = joinpath(data_folder, string(
        xlsx_name, ".", sheet_name, ".csv"
    ))
    
    tbl = XLSX.gettable(sheet, sheet_range;
        header=false,
        first_row=1, 
        infer_eltypes=false, 
        keep_empty_rows=true
    )
    global df  = DataFrame(tbl)                            # <- works via Tables.jl

    # first(df, 5)
    CSV.write(csv_file, df)

end


## ---.... .. - . .-  . -. . . - .  . . - -
## ---.... .. - . .-  . -. . . - .  . . - -
s
# Convert a Matrix (from XLSX.readdata) to a DataFrame with light header handling
function matrix_to_df(mat::AbstractMatrix)
    # Replace `nothing` with `missing` so CSV writes cleanly
    mat2 = map(x -> x === nothing ? missing : x, mat)
    if size(mat2, 1) == 0
        return DataFrame()
    end
    # If first row is all strings/symbols, use as header; otherwise auto names.
    firstrow = vec(mat2[1, :])
    headerish = all(x -> x isa AbstractString || x isa Symbol, firstrow)
    if headerish
        names = Symbol.(string.(firstrow))
        data  = size(mat2,1) > 1 ? mat2[2:end, :] : Array{Any}(undef, 0, size(mat2,2))
        return DataFrame(data, names)
    else
        names = Symbol.("Col", 1:size(mat2,2))
        return DataFrame(mat2, names)
    end
end

function read_sheet_to_df(xf::XLSX.XLSXFile, sheetname::String; range::Union{Nothing,String}=nothing)
    sh = xf[sheetname]
    if range === nothing
        # Easiest path: use Tables interface, which reads the used range and takes the
        # first row as header when appropriate (common case).
        try
            tbl = XLSX.readtable(sh)  # uses the sheet's used range
            return DataFrame(tbl)
        catch
            # Fallback: read the whole used range as raw matrix, then infer header
            ref = XLSX.get_ref(sh)  # e.g. "A1:D20"
            mat = XLSX.readdata(sh, ref)
            return matrix_to_df(mat)
        end
    else
        mat = XLSX.readdata(sh, range)
        return matrix_to_df(mat)
    end
end

function main()
    xlsx_path, opts = parse_args(ARGS)
    range = opts["range"] === nothing ? nothing : String(opts["range"])
    using_all = opts["all"] === true

    # Output path logic
    explicit_out = opts["out"] === nothing ? nothing : String(opts["out"])

    XLSX.openxlsx(xlsx_path) do xf
        sheets = XLSX.sheetnames(xf)
        # Resolve chosen sheets
        chosen_sheets::Vector{String}
        if using_all
            explicit_out !== nothing && error("--out is only for single-sheet mode.")
            chosen_sheets = sheets
        else
            sh_opt = opts["sheet"]
            if sh_opt === nothing
                chosen_sheets = [sheets[1]]
            else
                s = String(sh_opt)
                if occursin(r"^\d+$", s)
                    idx = parse(Int, s)
                    1 <= idx <= length(sheets) || error("Sheet index $idx out of bounds (1..$(length(sheets))).")
                    chosen_sheets = [sheets[idx]]
                else
                    s in sheets || error("Sheet \"$s\" not found. Available: $(join(sheets, \", \")).")
                    chosen_sheets = [s]
                end
            end
        end

        # Base name for auto output files
        base = replace(basename(xlsx_path), r"\.xlsx$"i => "")
        for (k, sname) in pairs(chosen_sheets)
            df = read_sheet_to_df(xf, sname; range=range)

            outpath = if explicit_out !== nothing
                explicit_out
            else
                suffix = using_all || length(chosen_sheets) > 1 || (opts["sheet"] !== nothing) ? "_$(sname)" : ""
                "$(base)$(suffix).csv"
            end

            CSV.write(outpath, df)
            println("✓ Wrote $(nrow(df)) rows × $(ncol(df)) cols to \"$outpath\" from sheet \"$sname\"" *
                    (range === nothing ? "" : " (range $range)"))
        end
    end
end

isinteractive() || main()
