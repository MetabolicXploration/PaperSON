src = """
x, y
0.28894016528253935, 0.028946967235725155
21.241245516372725, 0.02933170732145171
35.70409021148209, 0.015112953387066242
49.862317670638284, 0.01411444326281324
53.28479948638171, 1.0856390141186478
71.21749584990887, 2.866732393097978
92.70936533660233, 3.8633626576393842
99.82271165315306, 4.77776738975347
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





















