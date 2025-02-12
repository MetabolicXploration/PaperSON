src = """
x, y
6.251476571375545, 28.17022921625693
8.217063977306845, 63.55695768698362
8.613311618339711, 81.80355152688131
9.396944278111508, 111.40218491100163
9.861589778262973, 148.27991427343724
10.482070841132115, 156.55574557690417
10.862349959940005, 190.55963187235176
"""

# extract
let
    # clear lines
    lines = split(src, "\n"; keepempty = false)
    filter!(lines) do line
        contains(line, "x, y") && return false
        contains(line, ",") || return false
        return true
    end

    xs = Float64[]
    ys = Float64[]
    for line in lines
        xstr, ystr = split(line, ",")
        x = parse(Float64, xstr)
        y = parse(Float64, ystr)
        push!(xs, x)
        push!(ys, y)
    end
    # comma_lines = String[]

    # Copy
    clipboard(string(
        join(repr.(xs), ", "), 
        "\n\n", 
        join(repr.(ys), ", "), 
    ))
end
