src = """
x, y
3.622887848982237, 0.9869606205727337
4.064403993373234, 1.5813720125043336
4.645680938576918, 1.481531914034092
5.05055857783066, 1.4659201765056917
5.26359813188636, 1.458489416299088
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
